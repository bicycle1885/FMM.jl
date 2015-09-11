# concatenated chromosomes
immutable Genome
    seq::NucleotideSequence
    names::Vector{ASCIIString}
    offsets::Vector{Int}
end

function load_genome(fasta)
    seq = DNASequence()
    names = ASCIIString[]
    offsets = Int[]
    for chr in read(fasta, FASTA)
        push!(names, chr.name)
        push!(offsets, length(seq) + 1)
        seq *= chr.seq
    end
    return Genome(seq, names, offsets)
end

length(genome::Genome) = length(genome.seq)

@inline function getindex(genome::Genome, i::Integer)
    return genome.seq[i]
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


# precomputed SA ranges for each DNA kmer
immutable SATable{T,k}
    data::Vector{UnitRange{T}}
end

function getindex{T,k}(tbl::SATable{T,k}, kmer::Kmer{DNANucleotide,k})
    return tbl.data[convert(UInt64, kmer)+1]
end


immutable GenomeIndex{T,k}
    genome::Genome
    fmindex::FMIndex{2,T}
    sa_table::SATable{T,k}
end

# fasta: stream of chromosomes
# k: size of precomputed k-mer table
# r: sampling interval of SA values
function build_index(fasta, k::Int=12, r::Int=4)
    genome = load_genome(fasta)
    ivec = two_bits_encode(genome.seq)
    σ = 4  # A/C/G/T
    fmindex = FMIndex(ivec, σ, r=r)
    sa_table = build_sa_table(fmindex, Kmer{DNANucleotide,k})
    return GenomeIndex(genome, fmindex, sa_table)
end

function two_bits_encode(seq)
    n = length(seq)
    ivec = IntVector{2,UInt8}(n)
    for i in 1:n
        nuc = seq[i]
        if nuc == DNA_N
            # 'N's are replaced with 'A' (0x00)
            # to encode nucleotides with 2 bits
            ivec[i] = convert(UInt8, DNA_A)
        else
            ivec[i] = convert(UInt8, nuc)
        end
    end
    return ivec
end

function build_sa_table{n,k,T}(fmindex::FMIndex{n,T}, ::Type{Kmer{DNANucleotide,k}})
    # TODO there should be a much faster way using a stack
    len = 4^k
    sa_table = Vector{UnitRange{T}}(len)
    for i in 1:len
        kmer = convert(Kmer{DNANucleotide,k}, UInt64(i-1))
        sa_table[i] = FMIndices.sa_range(kmer, fmindex)
    end
    return SATable{T,k}(sa_table)
end
