# concatenated chromosomes
immutable Genome
    seq::DNASequence
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
    return Genome(seq, names, offsets)
end

Base.length(genome::Genome) = length(genome.seq)
Base.endof(genome::Genome)  = length(genome.seq)

@inline function Base.getindex(genome::Genome, i::Integer)
    return genome.seq[i]
end

function copy_nucleotides(genome::Genome, startpos::Integer, stoppos::Integer)
    len = abs(stoppos - startpos) + 1
    nucs = Vector{DNANucleotide}(len)
    step = startpos < stoppos ? 1 : -1
    j = startpos
    for i in 1:len
        nucs[i] = genome.seq[j]
        j += step
    end
    return nucs
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