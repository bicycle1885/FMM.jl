# mock
type Alignment
    score::Score
end

function call(::Type{Alignment})
    return Alignment(typemin(Score))
end

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
    alignment::Alignment
    # alignment score cache
    scores::Vector{Score}
    # unpacked read sequence cache
    rseq::Vector{DNANucleotide}
    # unpacked genome sequence cache
    gseqs::Vector{Vector{DNANucleotide}}
    # DNASeqs
    dnaseqs::Vector{DNASeq}
    function ReadState()
        new(
            DNASequence(), DNASequence(), [], [], 0, 0,
            SmallPQueue{Score,SeedHitExt}(5),
            false,
            Alignment(), [], [], [], [],
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
    rs.alignment = Alignment()
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
    return rs.alignment
end

function resize_scores!(rs::ReadState, size)
    resize!(rs.scores, size)
end

function reset_scores!(rs::ReadState)
    fill!(rs.scores, Score(0))
end

function prepare_alignment!(rs::ReadState, n)
    resize!(rs.gseqs, n)
    resize!(rs.dnaseqs, n)
    fill!(rs.gseqs, DNANucleotide[])
end


# seedhit iterator

immutable SeedHitIterator
    rs::ReadState
    forward::Bool
end

function each_forward_seedhit(rs::ReadState)
    return SeedHitIterator(rs, true)
end

function each_reverse_seedhit(rs::ReadState)
    return SeedHitIterator(rs, false)
end

Base.start(iter::SeedHitIterator) = 1
function Base.done(iter::SeedHitIterator, i)
    seedhits = iter.forward ? iter.rs.seedhits : iter.rs.seedhits′
    i > endof(seedhits)
end
function Base.next(iter::SeedHitIterator, i)
    seedhits = iter.forward ? iter.rs.seedhits : iter.rs.seedhits′
    seedhits[i], i + 1
end


# sequence unpacking

function unpack_left_read!(rs::ReadState, seedhit)
    read = isforward(seedhit) ? rs.read : rs.read′
    unpack_seq!(rs.rseq, read, seed_start(seedhit) - 1, 1)
end

function unpack_right_read!(rs::ReadState, seedhit)
    read = isforward(seedhit) ? rs.read : rs.read′
    unpack_seq!(rs.rseq, read, seed_stop(seedhit) + 1, endof(read))
end

function unpack_genome!(rs::ReadState, i, genome, startpos, stoppos)
    unpack_seq!(rs.gseqs[i], genome, startpos, stoppos)
    rs.dnaseqs[i] = DNASeq(rs.gseqs[i])
end

function unpack_seq!(dst, src, startpos, stoppos)
    len = abs(startpos - stoppos) + 1
    resize!(dst, len)
    step = startpos ≤ stoppos ? 1 : -1
    i = 0
    j = startpos
    @inbounds while i < len
        dst[i+=1] = src[j]
        j += step
    end
    return dst
end
