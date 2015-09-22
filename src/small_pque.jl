# small priority-queue-like data structure
type SmallPQueue{S,T}
    size::Int
    max_size::Int
    scores::Vector{S}
    elms::Vector{T}
    function SmallPQueue(max_size)
        q = new(0, max_size, [], [])
        sizehint!(q.scores, max_size)
        sizehint!(q.elms, max_size)
        return q
    end
end

Base.length(q::SmallPQueue) = q.size

isfull(q::SmallPQueue) = q.size == q.max_size

function enqueue!{S,T}(q::SmallPQueue{S,T}, elm::T, score::S)
    if isfull(q)
        for i in 1:q.size
            if score > q.scores[i]
                push_out!(q, i)
                q.scores[i] = score
                q.elms[i] = elm
                break
            end
        end
    else
        if q.size > 0 && score < q.scores[end]
            loc = q.size + 1
        else
            loc = 1
            for i in 1:q.size
                if score > q.scores[i]
                    loc = i
                    break
                end
            end
        end
        insert!(q.scores, loc, score)
        insert!(q.elms, loc, elm)
        q.size += 1
    end
    return q
end

function push_out!(q, from)
    for i in from:q.size-1
        q.scores[i+1] = q.scores[i]
        q.elms[i+1] = q.elms[i]
    end
end

function Base.empty!(q::SmallPQueue)
    q.size = 0
    empty!(q.scores)
    empty!(q.elms)
    return q
end

function Base.getindex(q::SmallPQueue, i::Integer)
    @assert 1 ≤ i ≤ q.size
    return q.scores[i], q.elms[i]
end
