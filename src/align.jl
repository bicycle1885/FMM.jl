using StatsBase: WeightVec, sample

# state and working space
type AlignmentState
    # read counter
    count::Int
    read_state::ReadState
end

# simple and fast consecutive range
# IntRange{true}  => forward
# IntRange{false} => backward
immutable IntRange{d}
    start::Int
    stop::Int
end

length(ir::IntRange{true})  = ir.stop - ir.start + 1
length(ir::IntRange{false}) = ir.start - ir.stop + 1
@inline getindex(ir::IntRange{true}, i::Integer)  = ir.start + (i - 1)
@inline getindex(ir::IntRange{false}, i::Integer) = ir.start - (i - 1)


function extend(read, s, genome, t, forward, current_best_score)
    @assert read[s] == genome[t]
    a = read
    b = genome
    gap_open   = Float32(3.0 + 5.0)
    gap_extend = Float32(3.0)
    if forward
        let
            ra = IntRange{true}(s, endof(a))
            rb = IntRange{true}(t, t + length(ra) - 1)
            return align(a, ra, b, rb, gap_open, gap_extend, current_best_score)
        end
    else
        let
            ra = IntRange{false}(s, 1)
            rb = IntRange{false}(t, t - (length(ra) - 1))
            return align(a, ra, b, rb, gap_open, gap_extend, current_best_score)
        end
    end
end

"""
Align a short read `a[ra]` to a genome sequence `b[rb]` using affine penalty.
"""
function align(a, ra, b, rb, gap_open, gap_extend, current_best_score)
    m = length(ra)
    n = length(rb)
    model = MyScore
    # best scores for each prefix combinations: (a[ra[1:i]], b[rb[1:j]])
    H = Bio.Align.DPMatrix{Float32}(m, n)
    E = Bio.Align.DPMatrix{Float32}(m, n)
    F = Bio.Align.DPMatrix{Float32}(m, n)
    Bio.Align.fitsize!(H, m, n)
    Bio.Align.fitsize!(E, m, n)
    Bio.Align.fitsize!(F, m, n)
    # Gotoh's algorithm
    H[0,0] = 0
    # copy sequence `a` to an unpacked vector for fast access
    tmp_a = Vector{DNANucleotide}(m)
    for i in 1:m
        a_ra_i = tmp_a[i] = a[ra[i]]
        H[i,0] = -(gap_open + gap_extend * (i - 1))
        E[i,0] = typemin(Float32)
    end
    best_score = typemin(Float32)
    @inbounds for j in 1:n
        b_rb_j = b[rb[j]]
        H[0,j] = -(gap_open + gap_extend * (j - 1))
        F[0,j] = typemin(Float32)
        for i in 1:m
            E[i,j] = smax(
                E[i,j-1] - gap_extend,
                H[i,j-1] - gap_open
            )
            F[i,j] = smax(
                F[i-1,j] - gap_extend,
                H[i-1,j] - gap_open
            )
            H[i,j] = max(
                E[i,j],
                F[i,j],
                H[i-1,j-1] + model[tmp_a[i],b_rb_j]
            )
        end
        best_score = max(H[m,j], best_score)
        if best_score < current_best_score
            # cannot improve the current best score
            return best_score
        end
    end
    return best_score
end

# faster and simpler `max` function, but not safe
@inline function smax(x, y)
    ifelse(x > y, x, y)
end

function alignment_score(read, sp, genome, sp′, k, current_best_score)
    fscore = extend(read, sp + k - 1, genome, sp′ + k - 1, true, current_best_score)
    bscore = extend(read, sp, genome, sp′, false, current_best_score)
    # NOTE: maching score is assumed to be zero
    return fscore + bscore
end

function align_read!{T,k}(rs::ReadState, index::GenomeIndex{T,k}, profile)
    read = rs.read
    slen = profile.seed_length
    interval = profile.seed_interval
    n_seeds = fld(length(read) - slen, interval) + 1
    for i in 1:n_seeds
        seed_range = (i-1)*interval+1:(i-1)*interval+slen
        seed = read[seed_range]
        if hasn(seed)
            # seeds containing N cannot be used to exact maching
            continue
        else
            if length(seed) ≥ k
                kmer = convert(Kmer{DNANucleotide,k}, seed[end-k+1:end])
                init_range::UnitRange{Int} = index.sa_table[kmer]
                sa_range = FMIndices.sa_range(seed[1:end-k], index.fmindex, init_range)
                #@assert sa_range == FMIndices.sa_range(seed, index.fmindex)
            else
                sa_range = FMIndices.sa_range(seed, index.fmindex)
            end
            if !isempty(sa_range)
                push!(rs, SeedHit(seed_range, sa_range))
            end
        end
    end

    if !hashit(rs)
        # no clue
        rs.isaligned = false
        return rs
    end

    # find best alignment from matching seeds
    best_score = typemin(Score)
    if n_hits(rs) ≤ profile.max_trials_per_seedhit
        # try all seed hits
        for seedhit in rs
            alignment_scores!(rs, seedhit, index, profile.score_params, profile.max_trials_per_seedhit)
            best_score = max(maximum(rs.scores), best_score)
        end
    else
        n_seedhits = length(rs.seedhits)
        wv = Vector{Float32}(n_seedhits)
        for i in 1:endof(rs.seedhits)
            wv[i] = 1 / count(rs.seedhits[i])
        end
        for seedhit in sample(rs.seedhits, WeightVec(wv), n_seedhits, replace=false)
            alignment_scores!(rs, seedhit, index, profile.score_params, profile.max_trials_per_seedhit)
            best_score = max(maximum(rs.scores), best_score)
        end
    end

    rs.isaligned = true
    # TODO backtracking to calculate alignment
    rs.alignment = Alignment(best_score)
    return rs
end

# calculate alignment scores of all hit locations
function alignment_scores!(rs::ReadState, seedhit::SeedHit, index, params, max_trials_per_seedhit)
    # allocate space for scores
    n_hits = count(seedhit)
    if n_hits ≤ max_trials_per_seedhit
        idx = [1:n_hits;]
    else
        # sample seed hit locations
        idx = sample(1:n_hits, max_trials_per_seedhit, replace=false)
        n_hits = max_trials_per_seedhit
    end
    resize_scores!(rs, n_hits)
    reset_scores!(rs)
    prepare_alignment!(rs, N_PAR)

    if seed_start(seedhit) > 1
        # align left side of the read
        unpack_read!(rs, 1, seed_start(seedhit) - 1)
        i = 0
        while i < n_hits
            start_i = i
            n_filled = 0
            # fill as many genome sequences as possible (i.e. up to N_PAR)
            while i < n_hits && n_filled < N_PAR
                loc = FMIndices.sa_value(seedhit[idx[i+=1]], index.fmindex) + 1
                startpos = loc - 1
                stoppos = startpos - length(rs.rseq)
                unpack_genome!(rs, n_filled += 1, index.genome, startpos, stoppos)
            end
            alignment_scores!(rs.tmp_scores, params, rs.rseq, rs.dnaseqs, n_filled)
            for j in 1:n_filled
                rs.scores[start_i+j] += rs.tmp_scores[j]
            end
        end
    end

    # add scores of exact matching
    for i in 1:n_hits
        rs.scores[i] += params.matching_score
    end

    if seed_stop(seedhit) < length(rs.read)
        # align right size of the read
        unpack_read!(rs, seed_stop(seedhit) + 1, length(rs.read))
        i = 0
        while i < n_hits
            start_i = i
            n_filled = 0
            # fill as many genome sequences as possible (i.e. up to N_PAR)
            while i < n_hits && n_filled < N_PAR
                loc = FMIndices.sa_value(seedhit[idx[i+=1]], index.fmindex) + 1
                startpos = loc + seed_length(seedhit)
                stoppos = startpos + length(rs.rseq)
                unpack_genome!(rs, n_filled += 1, index.genome, startpos, stoppos)
            end
            alignment_scores!(rs.tmp_scores, params, rs.rseq, rs.dnaseqs, n_filled)
            for j in 1:n_filled
                rs.scores[start_i+j] += rs.tmp_scores[j]
            end
        end
    end
end
