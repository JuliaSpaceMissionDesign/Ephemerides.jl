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
