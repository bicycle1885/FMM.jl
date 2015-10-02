# precomputed SA ranges for each DNA kmer
immutable SATable{T,k}
    data::Vector{UnitRange{T}}
end

function getindex{T,k}(tbl::SATable{T,k}, kmer::Kmer{DNANucleotide,k})
    return tbl.data[convert(UInt64, kmer)+1]
end


# reference genome, FM-Index, and precomputed suffix array table
immutable GenomeIndex{T,k}
    genome::Genome
    fmindex::FMIndex{2,T}
    sa_table::SATable{T,k}
end

# k: size of precomputed k-mer table
# r: sampling interval of SA values
function build_index(fastafile, k::Int=12, r::Int=4)
    genome = load_genome(fastafile)
    ivec = two_bits_encode(genome.seq)
    σ = 4  # A/C/G/T
    Mbp = 1000^2
    if length(genome) < 500Mbp
        program = :SuffixArrays
        mmap = false
    else
        program = :psascan
        mmap = true
    end
    fmindex = FMIndex(ivec, σ, r=r, program=program, mmap=mmap)
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
        sa_table[i] = FMIndexes.sa_range(kmer, fmindex)
    end
    return SATable{T,k}(sa_table)
end
