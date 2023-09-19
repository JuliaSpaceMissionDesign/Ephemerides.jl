
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
- `nd` -- `Int32` number of double components in each array summary
- `ni` -- `Int32` number of integer components in each array summary
- `fwd` -- `Int32` record number of initial summary record
- `bwd` -- `Int32` record number of final summary record
- `ffa` -- `Int32` first free address of the file 
- `name` -- `String` internal name of the file
- `lend` -- `Bool` true if the file was generated in little endian 

### References 
- [DAF Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)

### See Also 
See also [`DAF`](@ref) and [`EphemerisProvider`](@ref)
"""
struct DAFHeader
    nd::Int       
    ni::Int       
    fwd::Int      
    bwd::Int      
    ffa::Int      
    name::String    
    lend::Bool      
end

"""
    DAFHeader(record::Vector{UInt8})

Parse the header (file record) of the DAF file, i.e., the first physical record in a DAF 
which contains global information about the file.

### References 
- [DAF Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)

### See Also 
See also [`DAF`](@ref), [`parse_daf_comment`](@ref) and [`parse_daf_summaries`](@ref)
"""
function DAFHeader(record::Vector{UInt8})

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
    initial_record(head::DAFHeader)

Return the record number of the initial summary record in the DAF 
"""
@inline initial_record(head::DAFHeader) = head.fwd

"""
    final_record(head::DAFHeader)

Return the record number of the final summary record in the DAF 
"""
@inline final_record(head::DAFHeader) = head.bwd

"""
    free_address(head::DAFHeader)

Return the first free address in the file, i.e., the address at which the first element of 
the next array is to be added.
"""
@inline free_address(head::DAFHeader) = head.ffa

"""
    endian(head::DAFHeader)

Return `true` if the DAF file is in little-endian.
"""
@inline endian(head::DAFHeader) = head.lend

"""
    filename(head::DAFHeader)

Return the internal description of the DAF.
"""
@inline filename(head::DAFHeader) = head.name

"""
    summary_size(head::DAFHeader)

Compute the size of a single summary record of a DAF file, in bytes.
"""
@inline summary_size(head::DAFHeader) = 8*(head.nd + (head.ni + 1) ÷ 2)


"""
    DAFSegmentDescriptor

A container object to store both SPK and PCK descriptors information.

### Fields 
- `segtype` -- `Int32` SPK/PCK segment type
- `tstart` -- `Float64` initial segment type, in TDB seconds since J2000.0
- `tend` -- `Float64` final segment type, in TDB seconds since J2000.0
- `tid` -- `Int32` target object NAIF ID
- `cid` -- `Int32` center object NAIF ID
- `axesid` -- `Int32` reference axes ID. Defaults to -1 for PCKs
- `iaa` -- `Int32` initial array address
- `faa` -- `Int32` final array address

### References 
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [PCK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html)
"""
struct DAFSegmentDescriptor
    segtype::Int       
    tstart::Float64      
    tend::Float64       
    tid::Int           
    cid::Int          
    axesid::Int       
    iaa::Int          
    faa::Int          
end

""" 
    DAFSegmentDescriptor(summary::Vector{UInt8}, head::DAFHeader, isspk::Bool)

Generate an SPK or PCK descriptor by parsing a DAF summary.
"""
function DAFSegmentDescriptor(summary::Vector{UInt8}, head::DAFHeader, isspk::Bool)
    if isspk
        parse_spk_segment_descriptor(summary, endian(head))
    else 
        parse_pck_segment_descriptor(summary, endian(head))
    end
end

"""
    segment_type(desc::DAFSegmentDescriptor)

Return the SPK/PCK segment type.
"""
@inline segment_type(desc::DAFSegmentDescriptor) = desc.segtype

"""
    initial_time(desc::DAFSegmentDescriptor)

Return the initial epoch of the interval for which ephemeris data are contained in the 
segment, in seconds since J2000.0
"""
@inline initial_time(desc::DAFSegmentDescriptor) = desc.tstart

"""
    final_time(desc::DAFSegmentDescriptor)

Return the final epoch of the interval for which ephemeris data are contained in the 
segment, in seconds since J2000.0
"""
@inline final_time(desc::DAFSegmentDescriptor) = desc.tend

"""
    center(desc::DAFSegmentDescriptor)

Return the NAIF integer code for the reference object or axes for SPK and PCK, respectively.
"""
@inline center(desc::DAFSegmentDescriptor) = desc.cid

"""
    target(desc::DAFSegmentDescriptor)

Return the NAIF integer code for the target object or axes for SPK and PCK, respectively.
"""
@inline target(desc::DAFSegmentDescriptor) = desc.tid

"""
    axes(desc::DAFSegmentDescriptor)

Return the NAIF integer code for the reference axes. It is valid only for SPK files and 
defaults to -1 for PCKs. 
"""
@inline axes(desc::DAFSegmentDescriptor) = desc.axesid

"""
    initial_address(desc::DAFSegmentDescriptor)

Return the initial address of the segment array in the DAF.
"""
@inline initial_address(desc::DAFSegmentDescriptor) = desc.iaa

"""
    final_address(desc::DAFSegmentDescriptor)

Return the final address of the segment array in teh DAF.
"""
@inline final_address(desc::DAFSegmentDescriptor) = desc.faa

""" 
    parse_spk_segment_descriptor(summary::Vector{UInt8}, lend::Bool)

Create a [`DAFSegmentDescriptor`](@ref) object by parsing a binary SPK segment descriptor.

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

Create a [`DAFSegmentDescriptor`](@ref) object by parsing a binary PCK segment descriptor. 
A default value of -1 is used to fill the reference frame field. The target and center
fields are used for the actual target and center axes.

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


"""
    DAF 

Container to hold the information of NAIF's Double precision Array File (DAF). 

### Fields 
- `filepath` -- `String` system filepath of the DAF 
- `array` -- `Vector{UInt8}` binary content of the DAF
- `header` -- `DAFHeader` file record of the DAF
- `comment` -- `String` text within the DAF comment area 
- `ftype` -- `Int` file type, equals 1 for SPK and 2 for PCK
- `desc` -- DAF PCK/SPK segment descriptors
- `seglist` -- `SPKSegmentList` list of the SPK/PCK segments within the DAF

### References 
- [DAF Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)

### See Also 
See also [`DAFHeader`](@ref), [`Ephemerides.SPKSegmentList`](@ref) and [`EphemerisProvider`](@ref)
"""
struct DAF 
    filepath::String
    array::Vector{UInt8}
    header::DAFHeader
    comment::String 
    ftype::Int
    desc::Vector{DAFSegmentDescriptor}
    seglist::SPKSegmentList
end

"""
    DAF(filename::String)

Parse a DAF file and retrieve its content. 

!!! note
    This function does not initialise the segment types. That operation can be done in a 
    second moment through the [`initialise_segments!`](@ref) function.

!!! warning 
    Only DAF files with identification word equal to "DAF/SPK" or "DAF/PCK" are 
    accepted, otherwise an error is thrown.

### See Also 
See also [`initialise_segments!`](@ref).
"""
function DAF(filename::String)

    # Retrieve file content using a memory map
    if !isfile(filename)
        systemerror("opening file: $(repr(filename))", Int32(2))
    end

    array = mmap(filename, Vector{UInt8});

    # Retrieve DAF identification word
    locidw = get_string(array, 0, 8)
    if locidw == "DAF/SPK" 
        ftype = 1 
    elseif locidw == "DAF/PCK"
        ftype = 2
    else 
        throw(
            jEph.EphemerisError(
                "Invalid ephemeris file format. Only DAF/SPK and DAF/PCK files are accepted."
            )
        )
    end

    # Retrieve DAF file record and comment section
    header = DAFHeader(array)
    comment = parse_daf_comment(array, header)

    # Retrieve the DAF summaries
    summaries = parse_daf_summaries(array, header)

    # Create the DAF Segment descriptors
    desclist = DAFSegmentDescriptor[]

    for k in reverse(eachindex(summaries))
        # Parse the DAF segment descriptor from the summary binary content 
        push!(desclist, DAFSegmentDescriptor(summaries[k], header, ftype == 1))
    end
    
    # Initialise a temporarily empty SPK segment list 
    seglist = SPKSegmentList()
    DAF(filename, array, header, comment, ftype, desclist, seglist)
    
end

function Base.show(io::IO, daf::DAF)
    println(io, "DAF/$(is_spk(daf) ? "SPK" : "PCK")")
    println(io, " filepath = \"$(daf.filepath)\"")
    println(io, " header = $(daf.header)")
end

"""
    get_comment(daf::DAF)

Return the comment written in the DAF comment section. 
"""
@inline comment(daf::DAF) = daf.comment

"""
    get_header(daf::DAF)

Return the [`DAFHeader`](@ref) header of the DAF.
"""
@inline header(daf::DAF) = daf.header

""" 
    get_array(daf::DAF) 

Return the byte content of the DAF file.
"""
@inline array(daf::DAF) = daf.array

"""
    get_descriptors(daf::DAF) 

Return the SPK/PCK segment descriptors contained in the DAF.
"""
@inline descriptors(daf::DAF) = daf.desc

"""
    get_segment_list(daf::DAF)

Return the [`Ephemerides.SPKSegmentList`](@ref) list of segments stored in the DAF.
"""
@inline segment_list(daf::DAF) = daf.seglist

"""
    filepath(daf::DAF)
    
Return the system path of the DAF.
"""
@inline filepath(daf::DAF) = daf.filepath

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
    initial_record(daf::DAF)

Return the record number of the initial summary record in the DAF.
"""
@inline initial_record(daf::DAF) = initial_record(header(daf))

"""
    final_record(daf::DAF)

Return the record number of the final summary record in the DAF.
"""
@inline final_record(daf::DAF) = final_record(header(daf))

"""
    free_address(daf::DAF)

Return the first free address in the file, i.e., the address at which the first element of 
the next array is to be added.
"""
@inline free_address(daf::DAF) = free_address(header(daf))

"""
    endian(daf::DAF)

Return `true` if the DAF is in little-endian.
"""
@inline endian(daf::DAF) = endian(header(daf))

""" 
    parse_daf_comment(array::Vector{UInt8}, header::DAFHeader)

Retrieve the comment section of a binary DAF.

### References 
- [DAF Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)

### See Also 
See also [`DAF`](@ref), [`DAFHeader`](@ref) and [`parse_daf_summaries`](@ref)
"""
function parse_daf_comment(array::Vector{UInt8}, header::DAFHeader)

    cmt = "" 
    @inbounds for idx = 2:initial_record(header)-1

        # Get record and look for EOT byte 
        record = get_record(array, idx)
        eot = findfirst(c->c==0x04, record)
    
        if isnothing(eot)
            if idx == initial_record(header)-1 
                @warn "Could not find the EOT byte in the DAF comment."
                # throw(ArgumentError("Could not find the EOT byte in the DAF comment."))
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
    parse_daf_summaries(array::Vector{UInt8}, head::DAFHeader)

Parse the DAF binary content and retrieve all the summary records.

### References 
- [DAF Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)

### See Also 
See also [`DAF`](@ref), [`DAFHeader`](@ref) and [`parse_daf_comment`](@ref)
"""
function parse_daf_summaries(array::Vector{UInt8}, head::DAFHeader)

    # This function neglects the summaries record names! 
    summaries = Vector{UInt8}[]
    nc = summary_size(head)

    # Keep parsing summaries until next != 0
    next = Int(initial_record(head)) 
    while next != 0 

        record = get_record(array, next)

        # Get the control items of this summary record
        next = Int(get_float(record, 0, endian(head)))
        nsum = Int(get_float(record, 16, endian(head)))

        # update summaries with all those found in record
        for j = 1:nsum 
            push!(summaries, record[25+(j-1)*nc:24+j*nc])
        end

    end

    return summaries

end


"""
    initialise_segments!(daf::DAF)

Fill the [`Ephemerides.SPKSegmentList`](@ref) by initialising the SPK/PCK segments associated to all 
the descriptors stores within the DAF.

### See Also 
See also [`DAF`](@ref) and [`create_spk_segment`](@ref)
"""
function initialise_segments!(daf::DAF)

    for desc in descriptors(daf)
        # Create the SPK segment type associated to the above descriptor
        add_segment!(segment_list(daf), create_spk_segment(daf, desc)) 
    end

    nothing
end

"""
    create_spk_segment(daf::DAF, desc::DAFSegmentDescriptor)

Initialise an SPK segment according to the segment type defined in the 
[`DAFSegmentDescriptor`](@ref) `desc`.
"""
function create_spk_segment(daf::DAF, desc::DAFSegmentDescriptor)
    
    if !(segment_type(desc) in keys(SPK_SEGMENTLIST_MAPPING))
        throw(jEph.EphemerisError(
            "unsupported SPK segment type $(segment_type(desc)) found in $(filepath(daf))."
        ))
    end 

    mapped_spktype = SPK_SEGMENTLIST_MAPPING[segment_type(desc)]
    if mapped_spktype == 1
        SPKSegmentType1(daf, desc)

    elseif mapped_spktype == 2
        SPKSegmentType2(daf, desc)

    elseif mapped_spktype == 3 
        SPKSegmentType8(daf, desc)

    elseif mapped_spktype == 4 
        SPKSegmentType9(daf, desc)

    elseif mapped_spktype == 5
        SPKSegmentType19(daf, desc)

    else
        SPKSegmentType20(daf, desc)
    end
    
end





