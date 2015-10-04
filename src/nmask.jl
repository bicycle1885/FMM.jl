# compressed mask of ambiguous nucleotides
immutable NMask
    blockmask::SucVector
    blocks::Vector{UInt64}
    len::Int
end

function NMask(bv::BitVector)
    n = length(bv.chunks)
    blockmask = BitVector()
    blocks = Vector{UInt64}()
    for i in 1:n
        chunk = bv.chunks[i]
        if chunk == 0
            # no N in this block
            push!(blockmask, false)
        else
            push!(blockmask, true)
            push!(blocks, chunk)
        end
    end
    return NMask(blockmask, blocks, length(bv))
end

@inline Base.length(nmask::NMask) = nmask.len

@inline function Base.getindex(nmask::NMask, i::Integer)
    blockid = Int(i - 1) >> 6 + 1
    #if !nmask.blockmask[blockid]
    if !IndexableBitVectors.unsafe_getindex(nmask.blockmask, blockid)
        return false
    end
    @inbounds block = nmask.blocks[rank1(nmask.blockmask, blockid)]
    bitid = Int(i - 1) & 0b111111 + 1
    return ((block >> (bitid - 1)) & 1) == 1
end
