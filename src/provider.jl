
"""
    EphemerisProvider(file::String)
    EphemerisProvider(files::Vector{String})

Create an `EphemerisProvider` instance by loading a single or multiple binary ephemeris 
kernel files specified by `files`. Currently, only NAIF Double precision Array File (DAF)
kernels (i.e., SPK and PCK) are accepted.

### Example 
```julia-repl 
julia> eph = EphemerisProvider("PATH_TO_KERNEL")
EphemerisProvider([...])

julia> eph = EphemerisProvider(["PATH_TO_KERNEL_1", "PATH_TO_KERNEL_2"])
EphemerisProvider([])
```
"""
struct EphemerisProvider <: jEph.AbstractEphemerisProvider
    files::Vector{DAF}
    descriptors::Vector{DAFSegmentDescriptor} 
    linktable::Dict{Int, Dict{Int, Vector{SPKLink}}}
end

EphemerisProvider(files::AbstractString) = EphemerisProvider([files])

function EphemerisProvider(files::Vector{<:AbstractString})

    # Initial parsing of each DAF file 
    ndafs = length(files)
    dafs = Vector{DAF}(undef, ndafs)

    seg_descr = DAFSegmentDescriptor[]
    
    # Temporary SPK links table 
    linktable = Dict{Int, Dict{Int, Vector{SPKLink}}}()

    @inbounds for did = ndafs:-1:1
        # Retrieve the header, comment, and content of each DAF 
        dafs[did] = DAF(files[did])

        # Retrieve all the summary records of this DAF  
        summaries = parse_daf_summaries(dafs[did])

        # Parse the descriptors from the last to first to account for 
        # segment priority ordering
        for k in eachindex(reverse(summaries))

            # Parse the segment descriptor
            desc = DAFSegmentDescriptor(dafs[did], summaries[k])
            push!(seg_descr, desc)

            # Create the SPK segment and add it to the DAF SPK list
            seg = create_spk_segment(dafs[did], desc)
            add_segment!(get_segment_list(dafs[did]), seg)

            # Add the spk links
            add_spklinks!(linktable, dafs[did], desc, seg, did)
        end
    end

    return EphemerisProvider(dafs, seg_descr, linktable) 

end


"""
    create_spk_segment(daf::DAF, desc::DAFSegmentDescriptor)

Initialise an SPK segment according to the segment type defined in the 
[`DAFSegmentDescriptor`](@ref) `desc`.
"""
function create_spk_segment(daf::DAF, desc::DAFSegmentDescriptor)
    if desc.segtype in (1,)
        SPKSegmentType1(daf, desc)
    elseif desc.segtype in (2, )
        SPKSegmentType2(daf, desc)
    end
end