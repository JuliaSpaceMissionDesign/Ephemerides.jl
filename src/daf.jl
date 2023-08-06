
"""
    FTPSTR 

Validation string that guarantees the integrity of a DAF file. 

### References 
- [DAF Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)
"""
const FTPSTR = "FTPSTR:\r:\n:\r\n:\r\x00:\x81:\x10\xce:ENDFTP";


"""
    DAFHeader 

The DAF header, or file record, is the first physical record in a DAF and stores general 
information about the content of the file. 

### Fields 
- `nd` -- `Int32` number of double components in each array summar
- `ni` -- `Int32` number of integer components in each array summa
- `fwd` -- `Int32` record number of initial summary record
- `bwd` -- `Int32` record number of final summary record
- `ffa` -- `Int32` first free address of the file 
- `name` -- `String` internal name of the file
-- `lend` -- `Bool` true if the file was generated in little endian 

### References 
- [DAF Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)

### See Also 
See also [`DAF`](@ref) and [`EphemerisProvider`](@ref)
"""
struct DAFHeader
    nd::Int32       
    ni::Int32       
    fwd::Int32      
    bwd::Int32      
    ffa::Int32      
    name::String    
    lend::Bool      
end

"""
    DAF 

### Fields 
- `filepath` -- `String`
- `array` -- `Vector{UInt8}`
- `header` -- `DAFHeader`
- `comment` -- `String`
- `ftype` -- `Int`
- `seglist` -- `SPKSegmentList`

### References 
- [DAF Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)

### See Also 
See also [`DAFHeader`](@ref), [`SPKSegmentList`](@ref) and [`EphemerisProvider`](@ref)
"""
struct DAF 
    filepath::String
    array::Vector{UInt8}
    header::DAFHeader
    comment::String 
    ftype::Int       # 1 for SPK, 2 for PCK 
    seglist::SPKSegmentList
end

"""
    DAF(filename::String)

Parse a DAF file and retrieve its content.

!!! warning 
    Only DAF files with identification word equal to "DAF/SPK" or "DAF/PCK" are 
    accepted, otherwise an error is thrown.
"""
function DAF(filename::String)

    # Retrieve file content using a memory map
    array = mmap(filename, Vector{UInt8});

    # Retrieve DAF identification word
    locidw = get_string(array, 0, 8)
    if locidw == "DAF/SPK" 
        ftype = 1 
    elseif locidw == "DAF/PCK"
        ftype = 2
    else 
        throw(
            EphemerisError(
                "Invalid ephemeris file format. Only DAF/SPK and DAF/PCK files are accepted."
            )
        )
    end

    # Retrieve DAF file record and comment section
    header = parse_daf_header(array)
    comment = parse_daf_comment(array, header)

    # Initialise segment types list 
    seglist = SPKSegmentList()

    DAF(filename, array, header, comment, ftype, seglist)
    
end

"""
    get_comment(daf::DAF)

Return the comment written in the DAF comment section. 
"""
@inline get_comment(daf::DAF) = daf.comment

"""
    get_header(daf::DAF)

Return the [`DAFHeader`](@ref) header of the DAF file.
"""
@inline get_header(daf::DAF) = daf.header

""" 
    get_array(daf::DAF) 

Return the byte content of the DAF file.
"""
@inline get_array(daf::DAF) = daf.array

"""
    get_segment_list(daf::DAF)

Return the [`SPKSegmentList`](@ref) list of segments stored in the DAF file.
"""
@inline get_segment_list(daf::DAF) = daf.seglist

"""
    filepath(daf::DAF)
    
Return the system path of the DAF file.
"""
@inline filepath(daf::DAF) = daf.filepath

"""
    summary_size(daf::DAF)

Compute the size of a single summary record of a DAF file, in bytes.
"""
@inline summary_size(daf::DAF) = 8*(daf.header.nd + (daf.header.ni + 1) ÷ 2)

"""
    is_spk(daf::DAF)

Return `true` if the DAF stores SPK data.
"""
@inline is_spk(daf::DAF) = daf.ftype == 1

"""
    is_pck(daf::DAF)

Return `true` if the DAF stores PCK data.
"""
@inline is_pck(daf::DAF) = daf.ftype == 2


"""
    parse_daf_header(record::Vector{UInt8})

Parse the header (file record) of the DAF file, i.e., the first physical record in a DAF 
which contains global information about the file.

### References 
- [DAF Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)

### See Also 
See also [`DAF`](@ref), [`parse_daf_comment`](@ref) and [`parse_daf_summaries`](@ref)
"""
function parse_daf_header(record::Vector{UInt8})

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

""" 
    parse_daf_comment(array::Vector{UInt8}, header::DAFHeader)

Retrieve the comment section of a binary DAF file.

### References 
- [DAF Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)

### See Also 
See also [`DAF`](@ref), [`parse_daf_header`](@ref) and [`parse_daf_summaries`](@ref)
"""
function parse_daf_comment(array::Vector{UInt8}, header::DAFHeader)

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

""" 
    parse_daf_summaries

### References 
- [DAF Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)

### See Also 
See also [`DAF`](@ref), [`parse_daf_header`](@ref) and [`parse_daf_comment`](@ref)
"""
function parse_daf_summaries(daf::DAF)

    # This function neglects the summaries record names! 
    summaries = Vector{UInt8}[]
    nc = summary_size(daf)

    # Keep parsing summaries until next != 0
    next = Int(daf.header.fwd) 
    while next != 0 

        record = get_record(daf.array, next)

        # Get the control items of this summary record
        next = Int(get_float(record, 0, daf.header.lend))
        nsum = Int(get_float(record, 16, daf.header.lend))

        # update summaries with all those found in record
        for j = 1:nsum 
            push!(summaries, record[25+(j-1)*nc:24+j*nc])
        end

    end

    return summaries

end


"""
    DAFSegmentDescriptor

### Fields 
- `segtype` -- `Int32`
- `tstart` -- `Float64`
- `tend` -- `Float64`
- `tid` -- `Int32`
- `cid` -- `Int32`
- `axesid` -- `Int32`
- `iaa` -- `Int32`
- `faa` -- `Int32`

### References 
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [PCK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html)
"""
struct DAFSegmentDescriptor

    segtype::Int32      # Segment Type 

    tstart::Float64     # Initial segment time 
    tend::Float64       # Final segment time

    tid::Int32          # Target object ID 
    cid::Int32          # Center object ID  (only for SPK)
    axesid::Int32       # Reference axes ID (only for SPK)

    iaa::Int32          # Initial array address
    faa::Int32          # Final array address
end


""" 
    DAFSegmentDescriptor(daf::DAF, summary::Vector{UInt8})

"""
function DAFSegmentDescriptor(daf::DAF, summary::Vector{UInt8})
    if daf.ftype == 1
        parse_spk_segment_descriptor(summary, daf.header.lend)
    else 
        parse_pck_segment_descriptor(summary, daf.header.lend)
    end
end


""" 
    parse_spk_segment_descriptor(summary::Vector{UInt8}, lend::Bool)

### References 
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
"""
function parse_spk_segment_descriptor(summary::Vector{UInt8}, lend::Bool)

    # Get initial and final epoch in seconds past J2000
    tstart = get_float(summary, 0, lend)
    tend = get_float(summary, 8, lend)

    # Get target, center and axes NAIF IDs
    tid = get_int(summary, 16, lend)
    cid = get_int(summary, 20, lend)
    aid = get_int(summary, 24, lend)

    # Get SPK segment type 
    segtype = get_int(summary, 28, lend)

    # Get initial and final array adddresses
    iaa = get_int(summary, 32, lend)
    faa = get_int(summary, 36, lend)

    DAFSegmentDescriptor(segtype, tstart, tend, tid, cid, aid, iaa, faa)

end

""" 
    parse_pck_segment_descriptor(summary::Vector{UInt8}, lend::Bool)

### References 
- [PCK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html)
"""
function parse_pck_segment_descriptor(summary::Vector{UInt8}, lend::Bool)

    # Get initial and final epoch in seconds past J2000
    tstart = get_float(summary, 0, lend)
    tend = get_float(summary, 8, lend)

    # Get target and reference axes NAIF IDs
    tid = get_int(summary, 16, lend)
    cid = get_int(summary, 20, lend)

    # Get PCK segment type 
    segtype = get_int(summary, 24, lend)

    # Get initial and final array adddresses
    iaa = get_int(summary, 28, lend)
    faa = get_int(summary, 32, lend)
    
    # For PCK segments, the reference axis defaults to -1, whereas the center 
    # is the set of axis wrt which these pair is defined
    DAFSegmentDescriptor(segtype, tstart, tend, tid, cid, -1, iaa, faa)

end



