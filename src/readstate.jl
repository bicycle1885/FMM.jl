# alignment state of a read
type ReadState
    # FASTA/FASTQ record
    record
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
    alignment::Nullable{AlignedRead}
    function ReadState()
        new(
            nothing,
            DNASequence(), DNASequence(), [], [], 0, 0,
            IntervalHeap{SeedHitExt}(),
            Nullable()
        )
    end
end

function setrecord!(rs::ReadState, record)
    rs.record = record
    rs.read = record.seq
    rs.read′ = reverse_complement(record.seq)
    empty!(rs.seedhits)
    empty!(rs.seedhits′)
    rs.n_hits = 0
    rs.n_hits′ = 0
    empty!(rs.seedhit_queue)
    rs.alignment = Nullable()
    return rs
end

function setalignment!(rs::ReadState, aln)
    rs.alignment = Nullable(aln)
    return rs
end

function record(rs::ReadState)
    return rs.record
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
    clusters, densities = find_clusters(rs, true, index.fmindex)
    clusters′, densities′ = find_clusters(rs, false, index.fmindex)
    append!(clusters, clusters′)
    append!(densities, densities′)
    ord = sortperm(densities, rev=true)
    #@show densities[ord]
    hits = Vector{Int}()
    for i in ord
        maxsize = 0
        maxj = 0
        for j in clusters[i]
            size = count(j > 0 ? rs.seedhits[j] : rs.seedhits′[-j])
            if size > maxsize
                maxsize = size
                maxj = j
            end
        end
        push!(hits, maxj)
    end
    return PermutedSeedHitIterator(rs, hits)
end

function find_clusters(rs, forward, fmindex)
    seedhits = forward ? rs.seedhits : rs.seedhits′
    inds = Vector{Int}()
    locs = Vector{Int}()
    for (i, seedhit) in enumerate(seedhits)
        seedlocs = Vector{Int}(count(seedhit))
        for j in 1:count(seedhit)
            loc = FMIndexes.sa_value(seedhit[j], fmindex) + 1
            push!(inds, forward ? i : -i)
            push!(locs, loc)
            seedlocs[j] = loc
        end
        seedhits[i] = attach(seedhit, seedlocs)
    end
    n_seedhits = length(seedhits)
    ord = sortperm(locs)
    # a cluster holds seedhit indices
    clusters = Vector{Int}[]
    densities = Int[]
    i = 1
    while i ≤ endof(ord)
        startpos = locs[ord[i]]
        cluster = [inds[ord[i]]]
        density = 1
        i += 1
        while i ≤ endof(ord) && (locs[ord[i]] - startpos) ≤ readlen(rs)
            push!(cluster, inds[ord[i]])
            density += 1
            i += 1
        end
        push!(clusters, cluster)
        push!(densities, density)
    end
    return clusters, densities
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
    ord = sortperm(locs)
    # starting positions of clusters
    clusters = Int[]
    i = 1
    while i ≤ endof(ord)
        startpos = locs[ord[i]]
        push!(clusters, startpos)
        i += 1
        while i ≤ endof(ord) && locs[ord[i]] - startpos ≤ readlen(rs)
            i += 1
        end
    end
    return clusters
end
