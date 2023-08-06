
struct SPKLink 
    desc::DAFSegmentDescriptor # The actual descriptor (with modified interval values!)
    did::Int # DAF File ID containing the data 
    fid::Int # Field number of in the SPKSegmentList associated to this type
    lid::Int # Index in the `fid`-th field of the SPKSegmentList from which to retrieve the SPK Segment
    fct::Int # 1 or -1 to check whether the link must be reversed! 
end


""" 
    reverse_link(link::SPKLink)

Reverse the sign of the SPKLink
"""
function reverse_link(link::SPKLink)
    SPKLink(link.desc, link.did, link.fid, link.lid, -link.fct)
end

# TODO: check whether table should become implementation of dict?

"""
    add_spklinks!(table::Dict, desc::DAFSegmentDescriptor, seg::AbstractSPKSegment, did::Int)
"""
function add_spklinks!(table::Dict, daf::DAF, desc::DAFSegmentDescriptor, 
            seg::AbstractSPKSegment, did::Int)

    # NB: the mapping is [FROM][TO] -> SegmentDescriptors

    # Retrieve field index, and inner list index
    fid = spk_field(seg)
    lid = length(getfield(daf.seglist, fid))

    # Create the forward and backward SPKLink if not already available 
    f_map = get!(table, desc.cid, Dict{Int, Vector{SPKLink}}())
    b_map = get!(table, desc.tid, Dict{Int, Vector{SPKLink}}())

    f_link = SPKLink(desc, did, fid, lid, 1)

    # Populate with both the forward and backward links, initialising the 
    # SPKLink key if a link between the two bodies was yet to be found
    push!(get!(f_map, desc.tid, SPKLink[]), f_link)
    push!(get!(b_map, desc.cid, SPKLink[]), reverse_link(f_link))

end