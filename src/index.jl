# precomputed SA ranges for each DNA kmer
immutable SATable{T,k}
    data::Vector{UnitRange{T}}
end

function getindex{T,k}(tbl::SATable{T,k}, kmer::Kmer{DNANucleotide,k})
    i = convert(UInt64, kmer)
    return tbl.data[i+1]
end

immutable GenomeIndex{T,k}
    genome::NucleotideSequence
    fmindex::FMIndex{2,T}
    sa_table::SATable{T,k}
end

function build_index(fasta)
    genome = mapreduce(x -> x.seq, (*), read(fasta, FASTA))
    ivec = int_vector(genome)
    σ = 4  # A/C/G/T
    r = 4
    fmindex = FMIndex(ivec, σ, r)
    sa_table = build_sa_table(fmindex, Kmer{DNANucleotide,12})
    return GenomeIndex(genome, fmindex, sa_table)
end

function int_vector(genome)
    n = length(genome)
    ivec = IntVector{2,UInt8}(n)
    for i in 1:n
        nuc = genome[i]
        if nuc == DNA_N
            ivec[i] = 0x00
        else
            ivec[i] = convert(UInt8, nuc)
        end
    end
    return ivec
end

function build_sa_table{n,k,T}(fmindex::FMIndex{n,T}, ::Type{Kmer{DNANucleotide,k}})
    # TODO there is much faster way
    len = 4^k
    sa_table = Vector{UnitRange{T}}(len)
    for i in 1:len
        kmer = convert(Kmer{DNANucleotide,k}, UInt64(i-1))
        sa_table[i] = FMIndices.sa_range(kmer, fmindex)
    end
    return SATable{T,k}(sa_table)
end
