
"""
    InterpCache{T}
"""
struct InterpCache{T}
    buff::Vector{DiffCache{Vector{T}, Vector{T}}}
end

"""
    InterpCache{T}(len::Int, buffsize::Int)

Create an `InterpCache` instance of type `T` with `len` buffers, each of size `buffsize`.
"""
function InterpCache{T}(len::Int, buffsize::Int) where T
    return InterpCache{T}([DiffCache(zeros(T, buffsize)) for _ in 1:len])
end

# Julia APIs 
Base.getindex(c::InterpCache, i) = c.buff[i]
Base.length(c::InterpCache) = length(c.buff)

"""
    get_buffer(c::InterpCache, idx::int, x::Number)

Return the `idx`-th buffer from the corresponding `DiffCache` depending on the type of `x`. 
"""
@inline get_buffer(c::InterpCache, idx::Int, x::Number) = get_tmp(c[idx], x)