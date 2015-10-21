using StatsBase: WeightVec, sample, sample!

function run_alignment(profile::AlignmentProfile, index, read_file)
    format = endswith(read_file, ".fa") ? FASTA :
             endswith(read_file, ".fq") ? FASTQ :
             error("unknown format")
    reads = open(read_file, format)
    readstate = ReadState()
    n_reads = 0
    info("aligning reads")
    t = @elapsed for rec in reads
        setread!(readstate, rec.seq)
        align_read!(readstate, index, profile)
        n_reads += 1
        println(rec.name)
        if isaligned(readstate)
            aln = alignment(readstate)
            chr, loc = locus(index.genome, aln[2].startpos)
            println(chr, ':', loc)
            println(aln)
        else
            println("not aligned")
        end
        println()
    end
    info("finished: ", t, " s")
    info(@sprintf("%.1f", n_reads / t), " reads/s")
end

function align_read!{T,k}(rs::ReadState, index::GenomeIndex{T,k}, profile)
    seedlen = profile.seed_length
    interval = profile.seed_interval

    # sarch seed hits
    search_seed!(rs, true,  index, seedlen, interval)
    search_seed!(rs, false, index, seedlen, interval)

    if !hashit(rs)
        # no clue
        return rs
    end


    # find best alignment from matching seeds
    max_trials_per_seedhit = profile.max_trials_per_seedhit
    max_seedcut_multiplier = profile.max_seedcut_multiplier
    for seedhit in each_forward_seedhit(rs)
        if count(seedhit) ≤ max_trials_per_seedhit * max_seedcut_multiplier
            score_seed!(rs, seedhit, index, profile.score_model, max_trials_per_seedhit)
        end
    end
    for seedhit in each_reverse_seedhit(rs)
        if count(seedhit) ≤ max_trials_per_seedhit * max_seedcut_multiplier
            score_seed!(rs, seedhit, index, profile.score_model, max_trials_per_seedhit)
        end
    end

    if !has_aligned_seed(rs)
        return rs
    end

    aln = align_hit(rs, index.genome, profile.score_model)
    #set_alignment!(rs, aln)
    return rs
end

function search_seed!{T,k}(rs::ReadState, forward, index::GenomeIndex{T,k}, seedlen, interval)
    read = forward ? forward_read(rs) : reverse_read(rs)
    for seed_range in SeedIterator(1:endof(read), seedlen, interval)
        seed = read[seed_range]
        if hasn(seed)
            # seeds containing N cannot be used to exact maching
            continue
        end
        if length(seed) ≥ k
            # use the SA range table of k-mers
            kmer::DNAKmer{k} = seed[end-k+1:end]
            init_range::UnitRange{Int} = index.sa_table[kmer]
            sa_range = FMIndexes.sa_range(seed[1:end-k], index.fmindex, init_range)
        else
            sa_range = FMIndexes.sa_range(seed, index.fmindex)
        end
        if !isempty(sa_range)
            push!(rs, SeedHit(seed_range, sa_range, forward))
        end
    end
end

# calculate alignment scores of all hit locations
function score_seed!(rs::ReadState, seedhit::SeedHit, index, model, max_trials_per_seedhit)
    indices = Vector{Int}()

    n_hits = count(seedhit)
    if n_hits ≤ max_trials_per_seedhit
        resize!(indices, n_hits)
        copy!(indices, 1:n_hits)
    else
        # sample seed hit locations
        n_hits = max_trials_per_seedhit
        resize!(indices, n_hits)
        sample!(1:n_hits, indices, replace=false)
    end

    #offset = 10
    offset = 0
    read = isforward(seedhit) ? forward_read(rs) : reverse_read(rs)

    locs = Vector{Int}(n_hits)
    for i in 1:n_hits
        locs[i] = FMIndexes.sa_value(seedhit[indices[i]], index.fmindex) + 1
    end

    refseqs = Vector{seq_t}(n_hits)

    # align left sequences
    for i in 1:n_hits
        len = seed_start(seedhit) - 1
        offset = locs[i] - 2
        refseqs[i] = subseq(index.genome, len, offset, true)
    end

    left_scores = Vector{Score}(n_hits)
    right_scores = Vector{Score}(n_hits)

    alns = paralign_score(
        model.submat,
        model.gap_open_penalty,
        model.gap_extend_penalty,
        read[1:seed_start(seedhit)-1],
        refseqs
    )
    for i in 1:n_hits
        left_scores[i] = alns[i].score
    end

    # align right sequences
    for i in 1:n_hits
        len = length(read) - seed_stop(seedhit)
        offset = locs[i] + seed_length(seedhit) - 1
        refseqs[i] = subseq(index.genome, len, offset, false)
    end

    alns = paralign_score(
        model.submat,
        model.gap_open_penalty,
        model.gap_extend_penalty,
        read[seed_stop(seedhit)+1:end],
        refseqs
    )
    for i in 1:n_hits
        right_scores[i] = alns[i].score
    end

    for i in 1:n_hits
        seedhitext = SeedHitExt(
            seedhit,
            locs[i],
            left_scores[i],
            0,  # TODO: fix
            right_scores[i]
        )
        push!(rs, seedhitext)
    end
end

# unpack `src[startpos:stoppos)` into `dst`
# note that `stoppos` is exclusive!
function unpack_seq!(dst, src, startpos, stoppos)
    len = abs(startpos - stoppos)
    resize!(dst, len)
    step = startpos < stoppos ? 1 : -1
    i = 0
    j = startpos
    @inbounds while i < len
        dst[i+=1] = src[j]
        j += step
    end
    return dst
end

function align_hit(rs::ReadState, genome::Genome, affinegap)
    # get the best alignment seed
    seedhit = best_aligned_seed(rs)
    score = total_score(seedhit)
    @show score, seedhit
end

# pairwise-alignment wrapper
function pairalign(rseq, gseq, affinegap)
    submat = affinegap.submat
    gap_open_penalty = affinegap.gap_open_penalty
    gap_extend_penalty = affinegap.gap_extend_penalty
    H, E, F = PairwiseAlignment.affinegap_global_align(rseq, gseq, submat, gap_open_penalty, gap_extend_penalty)
    m = length(rseq)
    max_score = typemin(Score)
    max_score_col = 0
    for j in 0:length(gseq)
        if H[m+1,j+1] > max_score
            max_score = H[m+1,j+1]
            max_score_col = j
        end
    end
    rseq′, gseq′ = PairwiseAlignment.affinegap_global_traceback(rseq, gseq, H, E, F, (m, max_score_col), submat, gap_open_penalty, gap_extend_penalty)
    return AlignmentResult(H[m+1,max_score_col+1], rseq′, gseq′)
end

function combine_gapped_sequences!(gseq, left_gseq, seedhit, right_gseq)
    append!(gseq, reversed_counts(left_gseq))
    push_chars!(gseq, seed_length(seedhit))
    append!(gseq, counts(right_gseq))
    return gseq
end
