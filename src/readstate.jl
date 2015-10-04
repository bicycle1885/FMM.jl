#          |<-- seed hit -->|
# ~~~~~~~~~^^^^^^^^^^^^^^^^^^~~~~~~~~~~~
#          ^                ^
#          seed_start       seed_stop
immutable SeedHit
    # read range
    rrange::UnitRange{Int}
    # suffix array range
    sarange::UnitRange{Int}
    # wheather forward or reverse complement
    forward::Bool
end

seed_start(seedhit::SeedHit) = start(seedhit.rrange)
seed_stop(seedhit::SeedHit) = last(seedhit.rrange)
seed_length(seedhit::SeedHit) = length(seedhit.rrange)
count(seedhit::SeedHit) = length(seedhit.sarange)
getindex(seedhit::SeedHit, i::Integer) = seedhit.sarange[i]
isforward(seedhit::SeedHit) = seedhit.forward


#              |<-left->|<-- seed hit -->|<-right->|
# read:        ~~~~~~~~~^^^^^^^^^^^^^^^^^^~~~~~~~~~~
# genome: --------------------------------------------------
#                       ^
#                       location
immutable SeedHitExt
    seedhit::SeedHit
    # genome location
    location::Int
    # alignment scores
    lscore::Score  # left
    hscore::Score  # seedhit
    rscore::Score  # right
end

total_score(x::SeedHitExt) = x.lscore + x.hscore + x.rscore

# delegete
seed_start(seedhit::SeedHitExt) = seed_start(seedhit.seedhit)
seed_stop(seedhit::SeedHitExt) = seed_stop(seedhit.seedhit)
seed_length(seedhit::SeedHitExt) = seed_length(seedhit.seedhit)
isforward(seedhit::SeedHitExt) = isforward(seedhit.seedhit)


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
    best::SmallPQueue{Score,SeedHitExt}
    # alignment status
    isaligned::Bool
    # best alignment
    alignment::Nullable{AlignmentResult}
    function ReadState()
        new(
            DNASequence(), DNASequence(), [], [], 0, 0,
            SmallPQueue{Score,SeedHitExt}(5),
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
    empty!(rs.best)
    rs.isaligned = false
    rs.alignment = Nullable()
    #rs.alignment = Alignment()
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

function n_total_hits(rs::ReadState)
    rs.n_hits + rs.n_hits′
end

function hashit(rs::ReadState)
    rs.n_hits + rs.n_hits′ > 0
end

function isaligned(rs::ReadState)
    rs.isaligned
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
