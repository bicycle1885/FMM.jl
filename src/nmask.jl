# compressed mask of ambiguous nucleotides
immutable NMask
    blockmask::SucVector
    blocks::Vector{UInt64}
    len::Int
end

function Base.convert(::Type{NMask}, bv::BitVector)
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

function Base.convert(::Type{BitVector}, nmask::NMask)
    bv = BitVector()
    for i in 1:length(nmask)
        push!(bv, nmask[i])
    end
    return bv
end

@inline Base.length(nmask::NMask) = nmask.len

@inline function divrem64(i)
    return i >> 6, i & 0b111111
end

@inline function block_bit(i)
    d, r = divrem64(Int(i - 1))
    return d + 1, r + 1
end

function block_start(blockid)
    return (blockid - 1) << 6 + 1
end

@inline function Base.getindex(nmask::NMask, i::Integer)
    blockid, bitid = block_bit(i)
    #if !nmask.blockmask[blockid]
    if !IndexableBitVectors.unsafe_getindex(nmask.blockmask, blockid)
        return false
    end
    block = nmask.blocks[rank1(nmask.blockmask, blockid)]
    return ((block >> (bitid - 1)) & 1) == 1
end

function findnext(nmask::NMask, i::Integer)
    if i > length(nmask)
        return 0
    end
    blockid, bitid = block_bit(i)
    if nmask.blockmask[blockid]
        # try to find in the current block
        block = nmask.blocks[rank1(nmask.blockmask, blockid)]
        d = findnext_in_block(block, bitid)
        if d > 0
            # found in the block
            return 64 * (blockid - 1) + d
        end
    end
    # search in the following blocks
    blockid = search1(nmask.blockmask, blockid + 1)
    if blockid == 0
        return 0
    end
    block = nmask.blocks[rank1(nmask.blockmask, blockid)]
    d = findnext_in_block(block, 1)
    @assert d > 0
    return 64 * (blockid - 1) + d
end

function findnext_in_block(block, bitid)
    block = block >> (bitid - 1)
    return block == 0 ? 0 : bitid + trailing_zeros(block)
end

function findprev(nmask::NMask, i::Integer)
    if i ≤ 0
        return 0
    end
    blockid, bitid = block_bit(i)
    if nmask.blockmask[blockid]
        # try to find in the current block
        block = nmask.blocks[rank1(nmask.blockmask, blockid)]
        d = findprev_in_block(block, bitid)
        if d > 0
            # found in the block
            return 64 * (blockid - 1) + d
        end
    end
    # search in the following blocks
    blockid = rsearch1(nmask.blockmask, blockid - 1)
    if blockid == 0
        return 0
    end
    block = nmask.blocks[rank1(nmask.blockmask, blockid)]
    d = findprev_in_block(block, 64)
    @assert d > 0
    return 64 * (blockid - 1) + d
end

function findprev_in_block(block, bitid)
    block = block << (64 - bitid)
    return block == 0 ? 0 : bitid - leading_zeros(block)
end

function hasn_within(nmask::NMask, r::UnitRange{Int})
    return findnext(nmask, first(r)) in r
end

function make_nbitmap(nmask::NMask, r::UnitRange{Int})
    ns = falses(length(r))
    i = findnext(nmask, first(r))
    while i in r
        ns[i-first(r)+1] = true
        i = findnext(nmask, i + 1)
    end
    return ns
end
