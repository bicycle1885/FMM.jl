# mock
type Alignment
    score::Float64
end

function call(::Type{Alignment})
    return Alignment(-Inf)
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
end

seed_start(seedhit::SeedHit) = start(seedhit.rrange)
seed_stop(seedhit::SeedHit) = last(seedhit.rrange)
seed_length(seedhit::SeedHit) = length(seedhit.rrange)
count(seedhit::SeedHit) = length(seedhit.sarange)
getindex(seedhit::SeedHit, i::Integer) = seedhit.sarange[i]


# alignment state of a read
type ReadState
    # read
    read::DNASequence
    # exact seed hits
    seedhits::Vector{SeedHit}
    # the number of hits
    n_hits::Int
    # alignment status
    isaligned::Bool
    # best alignment
    alignment::Alignment
    # alignment score cache
    scores::Vector{Score}
    # temporary score cache
    tmp_scores::Vector{Score}
    # unpacked read sequence cache
    rseq::Vector{DNANucleotide}
    # unpacked genome sequence cache
    gseqs::Vector{Vector{DNANucleotide}}
    # DNASeqs
    dnaseqs::Vector{DNASeq}
    function ReadState()
        new(DNASequence(), [], 0, false, Alignment(), [], [], [], [], [])
    end
end

function setread!(rs::ReadState, read)
    rs.read = read
    empty!(rs.seedhits)
    rs.n_hits = 0
    rs.isaligned = false
    rs.alignment = Alignment()
    return rs
end

function push!(rs::ReadState, hit::SeedHit)
    push!(rs.seedhits, hit)
    rs.n_hits += count(hit)
    return rs
end

function n_hits(rs::ReadState)
    rs.n_hits
end

function hashit(rs::ReadState)
    rs.n_hits > 0
end

# iterate over seed hits
start(rs::ReadState) = 1
done(rs::ReadState, i) = i > endof(rs.seedhits)
next(rs::ReadState, i) = rs.seedhits[i], i + 1

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
    resize!(rs.tmp_scores, n)
    resize!(rs.gseqs, n)
    resize!(rs.dnaseqs, n)
    fill!(rs.gseqs, DNANucleotide[])
end

function unpack_read!(rs::ReadState, startpos, stoppos)
    unpack_seq!(rs.rseq, rs.read, startpos, stoppos)
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

function choice(rs::ReadState)
    @assert hashit(rs)
    # TODO: sample without replacement
    r = rand()
    s = 0.0
    for i in 1:endof(rs.seedhits)
        s += length(rs.seedhits[i].sarange) / rs.n_hits
        if s ≥ r
            return i
        end
    end
    return endof(rs.seedhits)
end
