# alignment state of a read
type ReadState
    # read
    read::DNASequence
    read′::DNASequence
    # exact seed hits
    seedhits::Vector{SeedHit}
    seedhits′::Vector{SeedHit}
    # the number of hits
    n_hits::Int
    n_hits′::Int
    # best alignments
    seedhit_queue::IntervalHeap{SeedHitExt}
    # alignment status
    isaligned::Bool
    # best alignment
    alignment::Nullable{AlignmentResult}
    function ReadState()
        new(
            DNASequence(), DNASequence(), [], [], 0, 0,
            IntervalHeap{SeedHitExt}(),
            false,
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
    rs.isaligned = true
    return rs
end

# forward / reverse complement
forward_read(rs::ReadState) = rs.read
reverse_read(rs::ReadState) = rs.read′

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
