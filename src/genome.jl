# concatenated chromosomes
immutable Genome
    seq::GenomicSequence
    names::Vector{ASCIIString}
    offsets::Vector{Int}
end

function load_genome(fastafile)
    seq = DNASequence()
    names = ASCIIString[]
    offsets = Int[]
    for chr in open(fastafile, FASTA)
        push!(names, chr.name)
        push!(offsets, length(seq) + 1)
        seq *= chr.seq
    end
    #push!(offsets, length(seq) + 1)
    return Genome(seq, names, offsets)
end

Base.length(genome::Genome) = length(genome.seq)
Base.endof(genome::Genome)  = length(genome.seq)

@inline function Base.getindex(genome::Genome, i::Integer)
    return genome.seq[i]
end

function subseq(genome::Genome, len::Int, offset::Int, reversed::Bool)
    return subseq(genome.seq, len, offset, reversed)
end

function unpack_seq!(dst::Vector, genome::Genome, startpos, stoppos)
    # delegate
    return unpack_seq!(dst, genome.seq, startpos, stoppos)
end

function locus(genome::Genome, i::Integer)
    if !(1 ≤ i ≤ length(genome))
        throw(BoundsError(i))
    end
    # run binary search to locate a chromosome containing the i-th nucleotide
    offsets = genome.offsets
    lo = 1
    hi = length(offsets) + 1
    while hi - lo > 1
        m = div(lo + hi, 2)
        if i < offsets[m]
            hi = m
        else
            lo = m
        end
    end
    return genome.names[lo], i - offsets[lo] + 1
end

function eachname(genome::Genome)
    return genome.names
end

function Base.length(genome::Genome, name)
    i = findfirst(genome.names, name)
    if i == endof(genome.offsets)
        stop = length(genome.seq) + 1
    else
        stop = genome.offsets[i+1]
    end
    return stop - genome.offsets[i]
end
