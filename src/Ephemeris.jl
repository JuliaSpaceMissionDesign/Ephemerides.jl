
using Mmap

include("spk/spktypes.jl")

include("daf.jl")
include("utils.jl")

include("spk/spk1.jl")
include("spk/spk2.jl")

struct SPKLink 
    desc::DAFSegmentDescriptor # The actual descriptor (with modified interval values!)
    fid::Int # File containing the data 
    sid::Int # Index in the j-th field of the spk segment list of the spk type in desc
    fct::Int # 1 or -1 to check whether the link must be reversed! 
end

""" 
    reverse_link(link::SPKLink)

Reverse the sign of the SPKLink
"""
function reverse_link(link::SPKLink)
    SPKLink(link.desc, link.fid, link.sid, -link.fct)
end


struct EphemerisProvider
    files::Vector{DAF}
    descriptors::Vector{DAFSegmentDescriptor} # TODO: remove me add SPKLink dict
end

EphemerisProvider(files::String) = EphemerisProvider([files])

function EphemerisProvider(files::Vector{String})

    # Initial parsing of each DAF file 
    nfiles = length(files)
    dafs = Vector{DAF}(undef, nfiles)

    seg_descr = DAFSegmentDescriptor[] # TODO: remove me
    
    # Temporary SPK links table 
    linktable = Dict{Int, Dict{Int, Vector{SPKLink}}}()

    @inbounds for fid = nfiles:-1:1
        # Retrieve the header, comment, and content of each DAF 
        dafs[fid] = DAF(files[fid])

        # Retrieve all the summary records of this DAF  
        summaries = parse_daf_summaries(dafs[fid])

        # Parse the descriptors 
        for k in eachindex(summaries)

            # Parse the segment descriptor
            desc = DAFSegmentDescriptor(dafs[fid], summaries[k], fid)
            push!(seg_descr, desc) # TODO: remove me

            # Add the spk links
            add_spklinks!(linktable, desc, fid)
        end
    end

    return EphemerisProvider(dafs, seg_descr), linktable

end

"""
    add_spklinks!(table::Dict, desc::DAFSegmentDescriptor)
"""
function add_spklinks!(table::Dict, desc::DAFSegmentDescriptor, fid::Int)

    # NB: the mapping is [FROM][TO] -> SegmentDescriptors

    # Create the forward and backward SPKLink if not already available 
    f_map = get!(table, desc.cid, Dict{Int, Vector{SPKLink}}())
    b_map = get!(table, desc.tid, Dict{Int, Vector{SPKLink}}())

    f_link = SPKLink(desc, fid, 0, 1)

    # Populate with both the forward and backward links, initialising the 
    # SPKLink key if a link between the two bodies was yet to be found
    push!(get!(f_map, desc.tid, SPKLink[]), f_link)
    push!(get!(b_map, desc.cid, SPKLink[]), reverse_link(f_link))

end

function optimise_links(links::Vector{SPKLink})

    # Initialise with the highest priority segment
    t_start = [links[1].desc.tstart]
    t_end   = [links[1].desc.tend]

    # new_links = SPKLink[links[1]] # The fid and spk settings need to be changed! 
    
    for link in links[2:end]
        
        # Retrieve segment start and end times 
        ts, te = link.desc.tstart, link.desc.tend

        j = findfirst(x->x > ts, t_start)
        k = findfirst(x->x < te, t_end)

        is_lb = !isnothing(j)
        is_ub = !isnothing(k) 

        if is_lb && is_ub 
            # This segment expands both ends, all the intermediate segments 
            # should be removed! 
            t_start, t_end = [ts], [te]
            
        elseif is_lb 
            

        elseif is_ub 

        else 
            # Not certain that the segment should be discarded!

        end

    end

    # return new_links
    return t_start, t_end 
end


kernels = ["/home/michele/spice/kernels/spk/de440.bsp",
           "/home/michele/spice/kernels/spk/de430.bsp"]

# using BenchmarkTools
eph, linktable = EphemerisProvider(kernels);

optimise_links(linktable[399][3])

eph1 = EphemerisProvider(kernels[1]);
eph2 = EphemerisProvider(kernels[2]);

# @benchmark eph = EphemerisProvider($kernels)

# @noinline function test_find(table::Dict, cid, tid)
#     table[cid][tid]
# end