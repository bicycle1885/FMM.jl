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


immutable SimpleSubstMatrix{T} <: AbstractSubstitutionMatrix{T}
    matching_score::T
    mismatching_score::T
end

Base.getindex(subst_matrix::SimpleSubstMatrix, x, y) = ifelse(x == y, subst_matrix.matching_score, subst_matrix.mismatching_score)

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
            score_seed!(rs, seedhit, index, profile.score_params, max_trials_per_seedhit)
        end
    end
    for seedhit in each_reverse_seedhit(rs)
        if count(seedhit) ≤ max_trials_per_seedhit * max_seedcut_multiplier
            score_seed!(rs, seedhit, index, profile.score_params, max_trials_per_seedhit)
        end
    end

    if isempty(rs.best)
        return rs
    end

    submat = SimpleSubstMatrix(
        profile.score_params.matching_score,
        profile.score_params.mismatching_score
    )
    affinegap = AffineGapScoreModel{Score}(
        submat,
        profile.score_params.gap_open_penalty,
        profile.score_params.gap_ext_penalty
    )
    aln = align_hit(rs, index.genome, affinegap)
    set_alignment!(rs, aln)
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

# global variables for parallel alignment
const indices = Vector{Int}()
const locs = Vector{Int}(N_PAR)
const left_scores  = Vector{Score}(N_PAR)
const right_scores = Vector{Score}(N_PAR)
const readseq = Vector{DNANucleotide}()
const genome_seqs_left  = Vector{Vector{DNANucleotide}}(N_PAR)
const genome_seqs_right = Vector{Vector{DNANucleotide}}(N_PAR)
const dna_seqs_left = Vector{DNASeq}(N_PAR)
const dna_seqs_right = Vector{DNASeq}(N_PAR)
for j in 1:N_PAR
    genome_seqs_left[j] = []
    genome_seqs_right[j] = []
    dna_seqs_left[j] = DNASeq()
    dna_seqs_right[j] = DNASeq()
end

# calculate alignment scores of all hit locations
function score_seed!(rs::ReadState, seedhit::SeedHit, index, params, max_trials_per_seedhit)
    global indices, locs, left_scores, right_scores, readseq, genome_seqs_left, genome_seqs_right, dna_seqs_left, dna_seqs_right

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

    offset = 10
    read = isforward(seedhit) ? forward_read(rs) : reverse_read(rs)

    i = 0
    while i < n_hits
        n_filled = 0
        while i < n_hits && n_filled < N_PAR
            loc = FMIndexes.sa_value(seedhit[indices[i+=1]], index.fmindex) + 1
            locs[n_filled+=1] = loc
            # fill genome sequences - left
            startpos = loc - 1
            stoppos = startpos - (seed_start(seedhit) - 1 + offset)
            unpack_seq!(genome_seqs_left[n_filled], index.genome, startpos, stoppos)
            # fill genome sequences - right
            startpos = loc + seed_length(seedhit)
            stoppos = startpos + (length(read) - seed_stop(seedhit) + offset)
            unpack_seq!(genome_seqs_right[n_filled], index.genome, startpos, stoppos)
        end

        # calculate alignemnt scores - left
        fill!(left_scores, 0)
        unpack_seq!(readseq, read, seed_start(seedhit) - 1, 0)
        for j in 1:n_filled
            dna_seqs_left[j] = genome_seqs_left[j]
        end
        alignment_scores!(left_scores, params, readseq, dna_seqs_left, n_filled)

        # calculate alignemnt scores - right
        fill!(right_scores, 0)
        unpack_seq!(readseq, read, seed_stop(seedhit) + 1, endof(read) + 1)
        for j in 1:n_filled
            dna_seqs_right[j] = genome_seqs_right[j]
        end
        alignment_scores!(right_scores, params, readseq, dna_seqs_right, n_filled)

        for j in 1:n_filled
            seedhitext = SeedHitExt(
                seedhit,
                locs[j],
                left_scores[j],
                params.matching_score * seed_length(seedhit),
                right_scores[j]
            )
            enqueue!(rs.best, seedhitext, total_score(seedhitext))
        end
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
    score, seedhit = rs.best[1]
    loc = seedhit.location
    read = isforward(seedhit) ? forward_read(rs) : reverse_read(rs)
    offset = 10
    # left
    left_rseq = Vector{DNANucleotide}()
    left_gseq = Vector{DNANucleotide}()
    unpack_seq!(left_rseq, read, seed_start(seedhit) - 1, 0)
    startpos = loc - 1
    stoppos = startpos - (seed_start(seedhit) - 1 + offset)
    unpack_seq!(left_gseq, genome, startpos, stoppos)
    left_aln = pairalign(left_rseq, left_gseq, affinegap)
    # right
    right_rseq = Vector{DNANucleotide}()
    right_gseq = Vector{DNANucleotide}()
    unpack_seq!(right_rseq, read, seed_stop(seedhit) + 1, length(read) + 1)
    startpos = loc + seed_length(seedhit)
    stoppos = startpos + (length(read) - seed_stop(seedhit) + offset)
    unpack_seq!(right_gseq, genome, startpos, stoppos)
    right_aln = pairalign(right_rseq, right_gseq, affinegap)
    # create whole alignment
    read′ = GappedSequence(read, 1)
    combine_gapped_sequences!(read′, left_aln[1], seedhit, right_aln[1])
    genome′ = GappedSequence(genome.seq, loc - n_chars(left_aln[2]))
    combine_gapped_sequences!(genome′, left_aln[2], seedhit, right_aln[2])
    return AlignmentResult(score, read′, genome′)
end

# pairwise-alignment wrapper
function pairalign(rseq, gseq, affinegap)
    submat = affinegap.subst_matrix
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
