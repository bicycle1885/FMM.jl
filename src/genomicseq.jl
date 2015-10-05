# economy representation of a genome/chromosome sequence
# since `N` nucleotides occupies small portion of a reference sequence,
# these nucleotides can be compactly compressed while keeping fast access.
immutable GenomicSequence
    data::Vector{UInt64}
    nmask::NMask
    len::Int
end

function Base.convert(::Type{GenomicSequence}, seq::DNASequence)
    @assert seq.part[1] == 1
    return GenomicSequence(copy(seq.data), NMask(seq.ns), length(seq))
end

function Base.convert(::Type{DNASequence}, seq::GenomicSequence)
    return NucleotideSequence{DNANucleotide}(seq.data, seq.nmask, 1:length(seq), false, false)
end

Base.length(seq::GenomicSequence) = seq.len

@inline function divrem32(i)
    return i >> 5, i & 0b11111
end

@inline function getnuc(data::Vector{UInt64}, i::Integer)
    d, r = divrem32(Int(i) - 1)
    @inbounds return convert(DNANucleotide, convert(UInt8, (data[d+1] >> 2r) & 0b11))
end

@inline function Base.getindex(seq::GenomicSequence, i::Integer)
    return seq.nmask[i] ? DNA_N : getnuc(seq.data, i)
end

# unpack `src[startpos:stoppos)` into `dst`
# note that `stoppos` is exclusive!
function unpack_seq!(dst::Vector, src::GenomicSequence, startpos, stoppos)
    len = abs(startpos - stoppos)
    resize!(dst, len)
    step = startpos < stoppos ? 1 : -1
    # unpack nucleotides without caring about ambiguous nucleotides
    i = 0
    j = startpos
    while i < len
        @inbounds dst[i+=1] = getnuc(src.data, j)
        j += step
    end
    # fix the unpacked nucleotides if there are ambiguous nucleotides
    if step > 0
        j = startpos - 1
        while 0 < (j = findnext(src.nmask, j + 1)) < stoppos
            dst[j-startpos+1] = DNA_N
        end
    else
        j = startpos + 1
        while (j = findprev(src.nmask, j - 1)) > stoppos
            dst[startpos-j+1] = DNA_N
        end
    end
    return dst
end
