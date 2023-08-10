using Mmap
using StaticArrays, BenchmarkTools

struct NotImplementedError <: Exception 
    msg::String
end 

function Base.showerror(io::IO, err::NotImplementedError)
    print(io, "NotImplementedError: $(err.msg)")
end

# Validation string to check whether the file is intact!
const FTPSTR = "FTPSTR:\r:\n:\r\n:\r\x00:\x81:\x10\xce:ENDFTP";

# Each record in the DAF is of 1024 bytes 
const DAF_RECORD_LENGTH = 1024; 

@inline @views function get_record(array, index::Integer)
    @inbounds return array[1+DAF_RECORD_LENGTH*(index-1):DAF_RECORD_LENGTH*index]
end

# Read from array the string at address with given length
function get_string(array, address::Integer, bytes::Integer)
    # address is in 0-index notation!
    @inbounds rstrip(String(@view(array[address+1:address+bytes])))
end

# Return true if the array has a little-endian format! 
function is_little_endian(array)
    endian = get_string(array, 88, 8);
    if endian == "LTL-IEEE"
        return true
    elseif endian == "BIG-IEEE"
        return false
    else 
        throw(ErrorException("The endiannes could not be recognised!"))
    end
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

# assume this is a DAF/SPK!
struct DAFHeader
    nd::Int32
    ni::Int32
    fwd::Int32 
    bwd::Int32
    ffa::Int32 
    name::String
    lend::Bool
end

function read_daf_header(record)

    # Check FTP validation string 
    ftpval = get_string(record, 699, 28)
    if ftpval != FTPSTR
        throw(ArgumentError("The DAF FTP validation string is not valid.\n"))
    end 

    # Get spk endiannes 
    lend = is_little_endian(record)
    
    # Internal name/description of the file
    name = get_string(record, 16, 60)

    # Number of double and integer components
    nd = get_int(record, 8, lend)
    ni = get_int(record, 12, lend)

    # Record number of initial and final summary records
    fwd = get_int(record, 76, lend)
    bwd = get_int(record, 80, lend)

    # First free address in the file if you need to write something!
    ffa = get_int(record, 84, lend)
    
    DAFHeader(nd, ni, fwd, bwd, ffa, name, lend)

end

# Extract the comment from a binary spk file! 
function read_daf_comment(array, header::DAFHeader)

    cmt = "" 
    @inbounds for idx = 2:header.fwd-1
        # Get record and look for EOT byte 
        record = get_record(array, idx)
        eot = findfirst(c->c==0x04, record)
    
        if isnothing(eot)
            if idx == header.fwd-1 
                throw(ArgumentError("Could not find the EOT byte in the DAF comment."))
            end 
            
            # Remove all the null characters between two successive records
            idx = findlast(c->c!=0x00, record)
            crec = isnothing(idx) ? record[1:end] : record[1:idx]
        else 
            crec = record[1:eot]
        end
        # Replace null characters with new lines and concatenate
        cmt *= replace(String(crec), '\0'=>'\n')
    end

    return cmt

end 

# Read spice SPK file! 
function read_spk(array)
    header = read_daf_header(array)
end

# Read spice binary PCK file!
function read_pck(array)
    header = read_daf_header(array)
end

function load_daf(array)

    # Get identification word! 
    locidw = get_string(array, 0, 8)

    if locidw == "NAIF/SPK" 
        
    elseif locidw == "DAF/SPK" || locidw == "NAIF/DAF"
        read_spk(array)
    elseif locidw == "DAF/PCK"
        read_pck(array)
    else
        throw(NotImplementedError("Unsupported SPICE kernel ($locidw).\n")) 
    end

end

# get control itmes of a record summary (without previous item address)
function get_summary_head(record, lend::Bool)
    next = Int(get_float(record, 0, lend))
    nsum = Int(get_float(record, 16, lend))
    return next, nsum
end


function get_daf_summaries(array, header::DAFHeader)
    # This function neglects the summaries record names! 

    summaries = Vector{UInt8}[]

    # Find summmary size!
    sumsize = header.nd + (header.ni + 1) ÷ 2
    nc = 8*sumsize 

    # Keep parsing summaries until next != 0
    next = Int(header.fwd) 
    while next != 0 

        record = get_record(array, next)
        next, nsum = get_summary_head(record, header.lend)

        # update summaries with all those found in record
        for j = 1:nsum 
            push!(summaries, record[25+(j-1)*nc:24+j*nc])
        end

    end

    return summaries

end

struct SPKSegmentHeader
    spktype::Int32 
    tstart::Float64 
    tend::Float64 
    tid::Int32 
    cid::Int32 
    axesid::Int32
    iaa::Int32 
    faa::Int32
end

# Parse each spk segment header!
function parse_spk_segment_header(summary, lend::Bool)

    # Get initial and final epoch in seconds past J2000
    tstart = get_float(summary, 0, lend)
    tend = get_float(summary, 8, lend)

    # Get target, center and axes NAIF IDs
    tid = get_int(summary, 16, lend)
    cid = get_int(summary, 20, lend)
    axesid = get_int(summary, 24, lend)

    # Get SPK segment type 
    spktype = get_int(summary, 28, lend)

    # Get initial and final array adddresses
    iaa = get_int(summary, 32, lend)
    faa = get_int(summary, 36, lend)

    SPKSegmentHeader(spktype, tstart, tend, tid, cid, axesid, iaa, faa)
end

