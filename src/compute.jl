
for (order, fun1, fun2) in zip(
    (1, 2, 3, 4),
    (:ephem_vector3, :ephem_vector6, :ephem_vector9, :ephem_vector12),
    (:spk_vector3, :spk_vector6, :spk_vector9, :spk_vector12)
)

    @eval begin 

        """
            $($fun1)(eph::EphemerisProvider, from::Int, to::Int, time::Number)

        Compute the $(3*$order)-elements state of one body (to) relative to another (from)
        at `time`, expressed in TDB/TCB seconds since J2000, in accordance with the kernel
        timescale.
        """
        function ($fun1)(eph::EphemerisProvider, from::Int, to::Int, time::Number)

            links = spk_links(eph)
            if haskey(links, from) && haskey(links[from], to)
                for link in links[from][to] 
                    if initial_time(link) <= time <= final_time(link)   
                        return factor(link)*$(fun1)(get_daf(eph, file_id(link)), link, time)
                    end
                end
            else 
                throw(
                    jEph.EphemerisError(
                        "ephemeris data for point with NAIFId $to with respect to point " * 
                        "$from is unavailable."
                    )
                )
            end
        
            throw(
                jEph.EphemerisError(
                    "ephemeris data for point with NAIFId $to with respect to point " *
                    "$from is not available at $(time) seconds since J2000."
                ),
            )

        end

        function ($fun1)(daf::DAF, link::SPKLink, time::Number)

            # Retrieve list and element link IDs
            lid = list_id(link)
            eid = element_id(link)
            
            if lid == 1 
                return $(fun2)(
                    daf, get_segment(segment_list(daf), 1, eid), time
                )
            elseif lid == 2
                return $(fun2)(
                    daf, get_segment(segment_list(daf), 2, eid), time
                )
            else
                return $(fun2)(
                    daf, get_segment(segment_list(daf), 3, eid), time
                )
            end
                
        end

    end

end