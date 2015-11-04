# alignment state of a read
type ReadState
    # read
    read::DNASequence
    read′::DNASequence
    # exact seed hits
    seedhits::Vector{SeedHit{Int}}
    seedhits′::Vector{SeedHit{Int}}
    # the number of hits
    n_hits::Int
    n_hits′::Int
    # best alignments
    seedhit_queue::IntervalHeap{SeedHitExt}
    # best alignment
    alignment::Nullable{AlignedSequence}
    function ReadState()
        new(
            DNASequence(), DNASequence(), [], [], 0, 0,
            IntervalHeap{SeedHitExt}(),
            Nullable()
        )
    end
end

function setread!(rs::ReadState, read)
    rs.read = read
    rs.read′ = reverse_complement(read)
    empty!(rs.seedhits)
    empty!(rs.seedhits′)
    rs.n_hits = 0
    rs.n_hits′ = 0
    empty!(rs.seedhit_queue)
    rs.alignment = Nullable()
    return rs
end

function set_alignment!(rs::ReadState, aln)
    rs.alignment = Nullable(aln)
    return rs
end

# forward / reverse complement
forward_read(rs::ReadState) = rs.read
reverse_read(rs::ReadState) = rs.read′

readlen(rs::ReadState) = length(rs.read)

function push!(rs::ReadState, hit::SeedHit)
    if isforward(hit)
        push!(rs.seedhits, hit)
        rs.n_hits += count(hit)
    else
        push!(rs.seedhits′, hit)
        rs.n_hits′ += count(hit)
    end
    return rs
end

function push!(rs::ReadState, seedhit::SeedHitExt)
    score = total_score(seedhit)
    push!(rs.seedhit_queue, seedhit)
    return rs
end

function best_aligned_seed(rs::ReadState)
    return maximum(rs.seedhit_queue)
end

function has_aligned_seed(rs::ReadState)
    return !isempty(rs.seedhit_queue)
end

function n_total_hits(rs::ReadState)
    rs.n_hits + rs.n_hits′
end

function hashit(rs::ReadState)
    rs.n_hits + rs.n_hits′ > 0
end

function isaligned(rs::ReadState)
    return !isnull(rs.alignment)
end

function alignment(rs::ReadState)
    if !isaligned(rs)
        error("this read is not aligned")
    end
    return get(rs.alignment)
end


# seedhit iterator

immutable SeedHitIterator
    rs::ReadState
    forward::Bool
end

each_forward_seedhit(rs::ReadState) = SeedHitIterator(rs, true)
each_reverse_seedhit(rs::ReadState) = SeedHitIterator(rs, false)

Base.start(iter::SeedHitIterator) = 1
function Base.done(iter::SeedHitIterator, i)
    seedhits = iter.forward ? iter.rs.seedhits : iter.rs.seedhits′
    i > endof(seedhits)
end
function Base.next(iter::SeedHitIterator, i)
    seedhits = iter.forward ? iter.rs.seedhits : iter.rs.seedhits′
    seedhits[i], i + 1
end

immutable PermutedSeedHitIterator
    rs::ReadState
    # positive: forward / negative: reverse
    ord::Vector{Int}
end

function each_permuted_seedhit(rs::ReadState)
    nforward = length(rs.seedhits)
    nreverse = length(rs.seedhits′)
    ord = Vector{Int}()
    for i in 1:nforward
        push!(ord, i)
    end
    for i in 1:nreverse
        push!(ord, -i)
    end
    shuffle!(ord)
    return PermutedSeedHitIterator(rs, ord)
end

Base.start(iter::PermutedSeedHitIterator) = 1
Base.done(iter::PermutedSeedHitIterator, i) = i > length(iter.ord)
function Base.next(iter::PermutedSeedHitIterator, i)
    idx = iter.ord[i]
    if idx > 0
        hit = iter.rs.seedhits[idx]
    else
        hit = iter.rs.seedhits′[-idx]
    end
    return hit, i + 1
end

function each_prioritized_seedhit(rs::ReadState, index)
    hits = Vector{Int}()
    append!(hits, 1:length(rs.seedhits))
    append!(hits, -1:-1:-length(rs.seedhits′))
    priority = Vector{Float64}()
    append!(priority, prioritize_seeds(rs, true,  index.fmindex))
    append!(priority, prioritize_seeds(rs, false, index.fmindex))
    ord = sortperm(priority, rev=true)
    return PermutedSeedHitIterator(rs, hits[ord])
end

function prioritize_seeds(rs, forward, fmindex)
    seedhits = forward ? rs.seedhits : rs.seedhits′
    inds = Vector{Int}()
    locs = Vector{Int}()
    for (i, seedhit) in enumerate(seedhits)
        seedlocs = Vector{Int}(count(seedhit))
        for j in 1:count(seedhit)
            loc = FMIndexes.sa_value(seedhit[j], fmindex) + 1
            push!(inds, i)
            push!(locs, loc)
            seedlocs[j] = loc
        end
        seedhits[i] = attach(seedhit, seedlocs)
    end
    n_seedhits = length(seedhits)
    #simmat = zeros(Int, n_seedhits, n_seedhits)
    priority = zeros(n_seedhits)
    ord = sortperm(locs)
    for i in 1:endof(ord)
        cent = locs[ord[i]]
        for j in i+1:endof(ord)
            if locs[ord[j]] - cent ≥ readlen(rs)
                break
            end
            priority[inds[ord[i]]] += 1
            priority[inds[ord[j]]] += 1
            #simmat[inds[ord[i]],inds[ord[j]]] += 1
            #simmat[inds[ord[j]],inds[ord[i]]] += 1
        end
    end
    #@show simmat
    return priority
end
