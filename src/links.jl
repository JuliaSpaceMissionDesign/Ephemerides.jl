
"""
    SPKLink

A link object to create a mapping between [`DAFSegmentDescriptor`](@ref) and its actual 
location within an [`EphemerisProvider`](@ref) object. 

### Fields 
- `desc` -- `DAFSegmentDescriptor` for the segment associated to this link
- `fid` -- `Int` index of the DAF containg the link data.
- `lid` -- `Int` field number in the [`SPKSegmentList`](@ref) for this segment type.
- `eid` -- `Int` index of the inner segment list that stores this SPK segment.
- `fct` -- `Int` 1 or -1 depending on whether the (from, to) directions must be reversed.

### See Also 
See also [`SPKLinkTable`](@ref), [`SPKSegmentList`](@ref) and [`add_spklinks!`](@ref).
"""
struct SPKLink 
    desc::DAFSegmentDescriptor 
    fid::Int 
    lid::Int 
    eid::Int 
    fct::Int 
end

"""
    descriptor(link::SPKLink)

Return the SPK/PCK segment descriptor associated to this link.
"""
@inline descriptor(link::SPKLink) = link.desc

"""
    file_id(link::SPKLink)

Return the DAF file index.
"""
@inline file_id(link::SPKLink) = link.fid

"""
    list_id(link::SPKLink)

Return the index of the list containing the segments of the given SPK/PCK type.
"""
@inline list_id(link::SPKLink) = link.lid

"""
    element_id(link::SPKLink)

Return the segment index in the inner SPK/PCK segment list.
"""
@inline element_id(link::SPKLink) = link.eid 

"""
    factor(link::SPKLink)

Return the direction multiplicative factor.
"""
@inline factor(link::SPKLink) = link.fct

"""
    initial_time(link::SPKLink)

Return the initial epoch of the interval for which ephemeris data are contained in the 
segment associated to this link, in seconds since J2000.0
"""
@inline initial_time(link::SPKLink) = initial_time(descriptor(link))

"""
    final_time(link::SPKLink)

Return the final epoch of the interval for which ephemeris data are contained in the 
segment associated to this link, in seconds since J2000.0
"""
@inline final_time(link::SPKLink) = final_time(descriptor(link))

""" 
    reverse_link(link::SPKLink)

Reverse the sign, i.e. change the sign of the multiplicative factor, of the link.
"""
function reverse_link(link::SPKLink)
    SPKLink(
        descriptor(link), file_id(link), list_id(link), element_id(link), -factor(link)
    )
end


"""
    SPKLinkTable 

Dictionary object providing all the [`SPKLink`](@ref) available between a set of (from, to)
objects
"""
SPKLinkTable = Dict{Int, Dict{Int, Vector{SPKLink}}} 

"""
    create_linktables(dafs::Vector{DAF})

Create the SPK and PCK [`SPKLinkTable`](@ref) for all the segments stores in the input DAFs.
"""
function create_linktables(dafs::Vector{DAF})

    spklinks, pcklinks = SPKLinkTable(), SPKLinkTable()
    for j in reverse(eachindex(dafs))
        add_spklinks!(is_spk(dafs[j]) ? spklinks : pcklinks, dafs[j], j)
    end

    return spklinks, pcklinks

end

"""
    add_spklinks!(table::SPKLinkTable, daf::DAF, fid::Int)

Insert in the input [`SPKLinkTable`](@ref) all the SPK or PCK links associated to 
the segment descriptors of the input DAF.
"""
function add_spklinks!(table::SPKLinkTable, daf::DAF, fid::Int)

    # Initialise the number of elements contained in each list
    nfields = fieldcount(SPKSegmentList)
    counter = zeros(nfields)

    for desc in descriptors(daf)

        segtype = segment_type(desc)

        # Get the field index of the list for this segment type
        lid = SPK_SEGMENTLIST_MAPPING[segtype]

        counter[lid] += 1

        # Create the forward and backward SPKLink if not already available 
        f_map = get!(table, center(desc), Dict{Int, Vector{SPKLink}}())
        b_map = get!(table, target(desc), Dict{Int, Vector{SPKLink}}())

        f_link = SPKLink(desc, fid, lid, counter[lid], 1)

        # Populate with both the forward and backward links, initialising the 
        # SPKLink key if a link between the two bodies was yet to be found
        push!(get!(f_map, target(desc), SPKLink[]), f_link)
        push!(get!(b_map, center(desc), SPKLink[]), reverse_link(f_link))
    end

end