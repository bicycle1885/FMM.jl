immutable SeedIterator
    start::Int
    readlen::Int
    seedlen::Int
    interval::Int
    function SeedIterator(range, seedlen, interval)
        @assert seedlen > 0 && interval > 0
        return new(first(range), length(range), seedlen, interval)
    end
end

Base.start(iter::SeedIterator) = 1
function Base.done(iter::SeedIterator, i)
    s = iter.start + (i - 1) * iter.interval
    return s + iter.seedlen - 1 > iter.readlen
end
function Base.next(iter::SeedIterator, i)
    s = iter.start + (i - 1) * iter.interval
    return s:(s+iter.seedlen-1), i + 1
end

function Base.length(iter::SeedIterator)
    return cld(iter.readlen - iter.seedlen + 1, iter.interval)
end
