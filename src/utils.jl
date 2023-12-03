"""
    DAF_RECORD_LENGTH

DAF record length, in bytes.

### References
- [DAF Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)
"""
const DAF_RECORD_LENGTH = 1024; 


"""
    get_record(array::Vector{UInt8}, index::Integer)

Retrieve a whole DAF record at position `index`.
"""
@inline @views function get_record(array, index::Integer)
    @inbounds return array[1+DAF_RECORD_LENGTH*(index-1):DAF_RECORD_LENGTH*index]
end


# Parsing Utils 
# =====================

"""
    is_little_endian(array::Vector{UInt8})

Return true if the array corresponds to the string indicating a little-endian format.
"""
function is_little_endian(array::Vector{UInt8})
    # TODO: add support for VAX-GFLT and VAX-DFLT
    endian = get_string(array, 88, 8);
    if endian == "LTL-IEEE"
        return true
    elseif endian == "BIG-IEEE"
        return false
    else 
        throw(ErrorException("The endiannes could not be recognised!"))
    end
end


# Read from array the string at address with given length
function get_string(array, address::Integer, bytes::Integer)
    # address is in 0-index notation!
    @inbounds rstrip(String(@view(array[address+1:address+bytes])))
end

@inline get_num(x::Number, lend::Bool) = lend ? htol(x) : hton(x)

function get_int(array, address::Integer, lend::Bool) 
    # address is in 0-index notation!
    ptr = unsafe_load(Ptr{Int32}(pointer(array, address+1)))
    get_num(ptr, lend)
end

function get_float(array, address::Integer, lend::Bool) 
    # address is in 0-index notation!
    ptr = unsafe_load(Ptr{Float64}(pointer(array, address+1)))
    get_num(ptr, lend)
end


# Vector Utils 
# =====================

vhat(u) = u/vnorm(u)
vnorm(u) = sqrt(vdot(u, u))
vsep(u, v) = acos(min(1, max(-1, vdot(vhat(u), vhat(v)))))

@inbounds vdot(u, v) = u[1]*v[1] + u[2]*v[2] + u[3]*v[3]
function vcross(u, v)
    return @inbounds SVector{3}(
        u[2]*v[3]-u[3]*v[2], u[3]*v[1]-u[1]*v[3], u[1]*v[2]-u[2]*v[1]
    )
end

function vrot(v1, k, θ)

    s, c = sincos(θ)
    u = vcross(k, v1)

    # Apply Rodrigues rotation formula to rotate `v1` around `k` of an angle `θ`
    return v1 + (1-c)*vcross(k, u) + s*u

end