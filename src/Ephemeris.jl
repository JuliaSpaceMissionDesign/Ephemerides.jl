
using Mmap
using StaticArrays

include("spk/spktypes.jl")

include("daf.jl")
include("utils.jl")

include("spk/spk1.jl")
include("spk/spk2.jl")

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


struct EphemerisProvider
    files::Vector{DAF}
    descriptors::Vector{DAFSegmentDescriptor} 
    linktable::Dict{Int, Dict{Int, Vector{SPKLink}}}
end

EphemerisProvider(files::String) = EphemerisProvider([files])

function EphemerisProvider(files::Vector{String})

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
            add_segment!(dafs[did].seglist, seg)

            # Add the spk links
            add_spklinks!(linktable, dafs[did], desc, seg, did)
        end
    end

    return EphemerisProvider(dafs, seg_descr, linktable) 

end

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

"""
    create_spk_segment(daf::DAF, desc::DAFSegmentDescriptor)
"""
function create_spk_segment(daf::DAF, desc::DAFSegmentDescriptor)
    if desc.segtype in (1,)
        SPKSegmentType1(daf, desc)
    elseif desc.segtype in (2, )
        SPKSegmentType2(daf, desc)
    end
end

function vector3(eph::EphemerisProvider, from::Int, to::Int, time::Number)

    # TODO: check that the path actually exists
    for link in eph.linktable[from][to] 
        if link.desc.tstart <= time <= link.desc.tend   
            return vector3(eph.files[link.fid], link, time)
        end
    end

    throw(ErrorException("Unable to find ephemeris data from $from to $to at time $time"));

end

function vector3(daf::DAF, link::SPKLink, time::Number)
    if link.desc.segtype in (2,)
        return link.fct*vector3(daf, daf.seglist.spk2[link.lid], link.desc, time)
    else 
        return link.fct*vector3(daf, daf.seglist.spk1[link.lid], link.desc, time)
    end
end



kernels = ["/home/michele/spice/kernels/spk/de440.bsp",
           "/home/michele/spice/kernels/spk/de430.bsp"]

# using BenchmarkTools
eph = EphemerisProvider(kernels);

eph1 = EphemerisProvider(kernels[1]);
eph2 = EphemerisProvider(kernels[2]);

vector3(eph, 399, 3, 0.0)

using BenchmarkTools
using FrameTransformations
@benchmark vector3($eph, 399, 3, 0.0)

using CALCEPH
using JSMDInterfaces.Ephemeris
using CalcephEphemeris

ephem = CalcephProvider(kernels[2])
prefetch(ephem)

compute(eph::Ephem,jd0::Float64,time::Float64,
   target::Integer,center::Integer,unit::Integer)

vector3(eph, 399, 3, 0.0)
compute(ephem, Float64(DJ2000), 0.0, 3, 399, useNaifId+unitKM+unitSec, 0)

y = @MVector zeros(3)
ephem_compute!(y, ephem, DJ2000, 0.0, 3, 399, 0)

jd0 = Float64(DJ2000)
unit = useNaifId+unitKM+unitSec

@benchmark vector3($eph, 399, 3, 0.0)
@benchmark ephem_compute!($y, $ephem, $DJ2000, 0.0, 3, 399, 0)


@code_warntype vector3(eph, 399, 3, 0.0)
@code_warntype vector4(eph, 399, 3, 0.0)

daf = eph.files[1];
link = eph.linktable[399][3][1]

@code_warntype vector3(daf, link, 0.0)
@code_warntype vector4(daf, link, 0.0)

@benchmark vector3($daf, $link, 0.0)
@benchmark vector4($daf, $link, 0.0)

# @benchmark eph = EphemerisProvider($kernels)

# @noinline function test_find(table::Dict, cid, tid)
#     table[cid][tid]
# end


# function optimise_links(links::Vector{SPKLink})

#     # Initialise with the highest priority segment
#     t_start = [links[1].desc.tstart]
#     t_end   = [links[1].desc.tend]

#     # new_links = SPKLink[links[1]] # The fid and spk settings need to be changed! 

#     for link in links[2:end]
        
#         # Retrieve segment start and end times 
#         ts, te = link.desc.tstart, link.desc.tend

#         j = findfirst(x->x > ts, t_start)
#         k = findfirst(x->x < te, t_end)

#         is_lb = !isnothing(j)
#         is_ub = !isnothing(k) 

#         println(is_lb, " ", is_ub)

#         if is_lb && is_ub 
#             # This segment expands both ends, all the intermediate segments 
#             # should be removed! 
#             t_start, t_end = [ts], [te]

#             # TODO: Here the new segment must be split into two parts, 
#             # from [ts, t_start], [t_end, te]
            
#         elseif is_lb 
#             println(j, " ", k)
#             # t_start = te < t_start[j+1] ? t_start[j+2:end] : t_start[j+1:end]
#             # t_end   = t_end[j+1:end]

#             # insert!(t_start, 1, ts)
#             # insert!(t_start, 1, te)

#         elseif is_ub 

#         else 
#             # Not certain that the segment should be discarded!

#         end

#     end

#     # return new_links
#     return t_start, t_end 
# end
