# TODO: check whether table should become implementation of dict?

"""
    SPKLink 
"""
struct SPKLink 
    desc::DAFSegmentDescriptor # The actual descriptor (with modified interval values!)
    did::Int # DAF File ID containing the data 
    fid::Int # Field number of in the SPKSegmentList associated to this type
    lid::Int # Index in the `fid`-th field of the SPKSegmentList from which to retrieve the SPK Segment
    fct::Int # 1 or -1 to check whether the link must be reversed! 
end

"""
    descriptor(link::SPKLink)

Return the [DAFSegmentDescriptor](@ref) associated to this link
"""
@inline descriptor(link::SPKLink) = link.desc

"""
    file_id(link::SPKLink)
"""
@inline file_id(link::SPKLink) = link.dif

"""
    field_id(link::SPKLink)
"""
@inline field_id(link::SPKLink) = link.fid

"""
    list_id(link::SPKLink)
"""
@inline list_id(link::SPKLink) = link.lid 

"""
    factor(link::SPKLink)

Return the direction multiplicative factor.
"""
@inline factor(link::SPKLink) = link.fct


""" 
    reverse_link(link::SPKLink)

Reverse the sign, i.e. change the sign of the multiplicative factor, of the link.
"""
function reverse_link(link::SPKLink)
    SPKLink(link.desc, link.did, link.fid, link.lid, -link.fct)
end

"""
    add_spklinks!(table::Dict, desc::DAFSegmentDescriptor, seg::AbstractSPKSegment, did::Int)
"""
function add_spklinks!(table::Dict, daf::DAF, desc::DAFSegmentDescriptor, 
            seg::AbstractSPKSegment, did::Int)

    # NB: the mapping is [FROM][TO] -> SegmentDescriptors

    # Retrieve field index, and inner list index
    fid = spk_field(seg)
    lid = length(getfield(get_segment_list(daf), fid))

    # Create the forward and backward SPKLink if not already available 
    f_map = get!(table, center(desc), Dict{Int, Vector{SPKLink}}())
    b_map = get!(table, target(desc), Dict{Int, Vector{SPKLink}}())

    f_link = SPKLink(desc, did, fid, lid, 1)

    # Populate with both the forward and backward links, initialising the 
    # SPKLink key if a link between the two bodies was yet to be found
    push!(get!(f_map, target(desc), SPKLink[]), f_link)
    push!(get!(b_map, center(desc), SPKLink[]), reverse_link(f_link))

end