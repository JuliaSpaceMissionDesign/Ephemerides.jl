
# Each record in the DAF is of 1024 bytes 
const DAF_RECORD_LENGTH = 1024; 

"""
    get_record(array::Vector{UInt8}, index::Integer)

Retrieve a whole DAF record at position `index`.
"""
@inline @views function get_record(array, index::Integer)
    @inbounds return array[1+DAF_RECORD_LENGTH*(index-1):DAF_RECORD_LENGTH*index]
end

"""
    is_little_endian(array::Vector{UInt8})

Return true if the array corresponds to the string indicating a little-endian format.
"""
function is_little_endian(array::Vector{UInt8})
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

# TODO: overload di get_flow in modo da non dover specificare ogni volta lend ma con DAF type