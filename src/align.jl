# alignment profile
immutable Profile
    seed_length::Int
    seed_interval::Int
    max_consective_fails::Int
    function Profile(;kwargs...)
        local seed_length, seed_interval, max_consective_fails
        for (kw, arg) in kwargs
            if kw === :seed_length
                seed_length = arg
            elseif kw === :seed_interval
                seed_interval = arg
            elseif kw === :max_consective_fails
                max_consective_fails = arg
            else
                error("unknown parameter: $kw")
            end
        end
        return new(seed_length, seed_interval, max_consective_fails)
    end
end

# mock
type Alignment
    score::Float64
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


# simple scoring model
immutable MyScoreModel; end
#@inline Base.getindex(::MyScoreModel, ::Character, ::Type{GAP}) = -3.0
#@inline Base.getindex(::MyScoreModel, ::Type{GAP}, ::Character) = -3.0
@inline Base.getindex(::MyScoreModel, x::Character, y::Character) = ifelse(x === y, 0.0, -3.0)
const MyScore = MyScoreModel()

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
            #return best_score
        end
    end
    return best_score
end

@inline function smax(x, y)
    ifelse(x > y, x, y)
end

function alignment_score(read, sp, genome, sp′, k, current_best_score)
    fscore = extend(read, sp + k - 1, genome, sp′ + k - 1, true, current_best_score)
    bscore = extend(read, sp, genome, sp′, false, current_best_score)
    # NOTE: maching score is assumed to be zero
    return fscore + bscore
end

function align_read{T,k}(read, index::GenomeIndex{T,k}, profile, ranges, weight)
    slen = profile.seed_length
    interval = profile.seed_interval
    n_seeds = fld(length(read) - slen, interval) + 1
    resize!(ranges, n_seeds)
    n_hits = 0
    for i in 1:n_seeds
        seed = read[(i-1)*interval+1:(i-1)*interval+slen]
        if hasn(seed)
            # seeds containing N cannot be used to exact maching
            ranges[i] = 0:-1
        else
            kmer = convert(Kmer{DNANucleotide,k}, seed[end-k+1:end])
            init_range = convert(UnitRange{Int}, index.sa_table[kmer])
            range = FMIndices.sa_range(seed[1:end-k], index.fmindex, init_range)
            #range′ = FMIndices.sa_range(seed, index.fmindex)
            #@assert range == range′
            ranges[i] = range
            n_hits += length(range)
        end
    end
    if n_hits == 0
        # no clue
        return Nullable{Alignment}()
    end
    best_score = -Inf
    if n_hits ≤ profile.max_consective_fails
        for i in 1:endof(ranges)
            range = ranges[i]
            for j in 1:length(range)
                sp = (i - 1) * interval + 1
                sp′ = FMIndices.sa_value(range[j], index.fmindex) + 1
                score = alignment_score(read, sp, index.genome, sp′, slen, best_score)
                best_score = max(score, best_score)
            end
        end
    else
        resize!(weight, n_seeds)
        for i in 1:n_seeds
            weight[i] = isempty(ranges[i]) ? 0.0 : length(ranges[i]) / n_hits
        end
        @assert sum(weight) ≈ 1.0
        seed_ext = profile.max_consective_fails
        while seed_ext > 0
            # TODO sampling without replacement
            # select seed
            i = choice(weight)
            # select seed hit
            range = ranges[i]
            j = rand(1:length(range))
            # run seed extension
            sp = (i - 1) * interval + 1
            sp′ = FMIndices.sa_value(range[j], index.fmindex) + 1
            score = alignment_score(read, sp, index.genome, sp′, slen, best_score)
            if score > best_score
                best_score = score
                seed_ext = profile.max_consective_fails
            else
                # seed extension failed
                seed_ext -= 1
            end
        end
    end
    return Nullable(Alignment(best_score))
end

function choice(weight)
    r = rand()
    s = 0.0
    for i in 1:endof(weight)
        s += weight[i]
        if s ≥ r
            return i
        end
    end
    return endof(weight)
end
