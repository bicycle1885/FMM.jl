using StatsBase: WeightVec, sample

function align_read!{T,k}(rs::ReadState, index::GenomeIndex{T,k}, profile)
    seedlen = profile.seed_length
    interval = profile.seed_interval

    # sarch seed hits
    search_seed!(rs, true,  index, seedlen, interval)
    search_seed!(rs, false, index, seedlen, interval)

    if !hashit(rs)
        # no clue
        rs.isaligned = false
        return rs
    end

    # find best alignment from matching seeds
    for seedhit in each_forward_seedhit(rs)
        score_seed!(rs, seedhit, index, profile.score_params, profile.max_trials_per_seedhit)
    end
    for seedhit in each_reverse_seedhit(rs)
        score_seed!(rs, seedhit, index, profile.score_params, profile.max_trials_per_seedhit)
    end

    # TODO backtracking to calculate alignment
    rs.isaligned = true
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
            kmer = convert(Kmer{DNANucleotide,k}, seed[end-k+1:end])
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
function score_seed!(rs::ReadState, seedhit::SeedHit, index, params, max_trials_per_seedhit)
    n_hits = count(seedhit)
    if n_hits ≤ max_trials_per_seedhit
        indices = [1:n_hits;]
    else
        # sample seed hit locations
        indices = sample(1:n_hits, max_trials_per_seedhit, replace=false)
        n_hits = max_trials_per_seedhit
    end

    locs = Vector{Int}(N_PAR)
    left_scores  = Vector{Score}(N_PAR)
    right_scores = Vector{Score}(N_PAR)
    genome_seqs_left  = Vector{Vector{DNANucleotide}}(N_PAR)
    genome_seqs_right = Vector{Vector{DNANucleotide}}(N_PAR)
    dna_seqs_left = Vector{DNASeq}(N_PAR)
    dna_seqs_right = Vector{DNASeq}(N_PAR)
    for j in 1:N_PAR
        genome_seqs_left[j] = []
        genome_seqs_right[j] = []
        dna_seqs_left[j] = DNASeq()
        dna_seqs_right[j] = DNASeq()
    end

    readlen_left = seed_start(seedhit) - 1
    readlen_right = length(rs.read) - seed_stop(seedhit)
    offset = 0

    i = 0
    while i < n_hits
        n_filled = 0
        while i < n_hits && n_filled < N_PAR
            loc = FMIndexes.sa_value(seedhit[indices[i+=1]], index.fmindex) + 1
            n_filled += 1
            locs[n_filled] = loc
            # fill genome sequences - left
            startpos = loc - 1
            stoppos = startpos - (readlen_left - 1) - offset
            unpack_seq!(genome_seqs_left[n_filled], index.genome, startpos, stoppos)
            # fill genome sequences - right
            startpos = loc + seed_length(seedhit)
            stoppos = startpos + (readlen_right - 1) + offset
            unpack_seq!(genome_seqs_right[n_filled], index.genome, startpos, stoppos)
        end
        # calculate alignemnt scores - left
        fill!(left_scores, 0)
        if seed_start(seedhit) > 1
            unpack_left_read!(rs, seedhit)
            for j in 1:n_filled
                dna_seqs_left[j] = genome_seqs_left[j]
            end
            alignment_scores!(left_scores, params, rs.rseq, dna_seqs_left, n_filled)
        end
        # calculate alignemnt scores - right
        fill!(right_scores, 0)
        if seed_stop(seedhit) < length(rs.read)
            unpack_right_read!(rs, seedhit)
            for j in 1:n_filled
                dna_seqs_right[j] = genome_seqs_right[j]
            end
            alignment_scores!(right_scores, params, rs.rseq, dna_seqs_right, n_filled)
        end
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
