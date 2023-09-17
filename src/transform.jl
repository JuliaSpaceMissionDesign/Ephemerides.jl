
export ephem_vector3, ephem_vector6, ephem_vector9, ephem_vector12, 
       ephem_rotation3, ephem_rotation6, ephem_rotation9, ephem_rotation12

for (order, pfun1, afun1, fun1) in zip(
    (1, 2, 3, 4),
    (:ephem_vector3, :ephem_vector6, :ephem_vector9, :ephem_vector12),
    (:ephem_rotation3, :ephem_rotation6, :ephem_rotation9, :ephem_rotation12),
    (:spk_vector3, :spk_vector6, :spk_vector9, :spk_vector12)
)

    @eval begin 

        """
            $($pfun1)(eph::EphemerisProvider, from::Int, to::Int, time::Number)

        Compute the $(3*$order)-elements state of one body (to) relative to another (from)
        at `time`, expressed in TDB/TCB seconds since J2000, in accordance with the kernel
        timescale.
        """
        function ($pfun1)(eph::EphemerisProvider, from::Int, to::Int, time::Number)

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

        """
            $($afun1)(eph::EphemerisProvider, from::Int, to::Int, time::Number)

        Compute the $(3*$order)-elements orientation angles of one set of axes (to) relative 
        to another (from) at `time`, expressed in TDB/TCB seconds since J2000, in accordance 
        with the kernel timescale.
        """
        function ($afun1)(eph::EphemerisProvider, from::Int, to::Int, time::Number)
            # TODO: check whether this -factor is also valid for orientation angles!
            links = pck_links(eph)
            if haskey(links, from) && haskey(links[from], to)
                for link in links[from][to] 
                    if initial_time(link) <= time <= final_time(link)   
                        return factor(link)*$(fun1)(get_daf(eph, file_id(link)), link, time)
                    end
                end
            else 
                throw(
                    jEph.EphemerisError(
                        "ephemeris data for axes with NAIFId $to with respect to axes " * 
                        "$from is unavailable."
                    )
                )
            end
        
            throw(
                jEph.EphemerisError(
                    "ephemeris data for axes with NAIFId $to with respect to axes " *
                    "$from is not available at $(time) seconds since J2000."
                ),
            )

        end

        function ($fun1)(daf::DAF, link::SPKLink, time::Number)

            # Retrieve list and element link IDs
            lid = list_id(link)
            eid = element_id(link)
            
            # Use binary search to reduce the time spent within the if\cycle
            if lid < 4  
                if lid == 1 
                    return $(fun1)(daf, get_segment(segment_list(daf), 1, eid), time)
                elseif lid == 2
                    return $(fun1)(daf, get_segment(segment_list(daf), 2, eid), time)
                else
                    return $(fun1)(daf, get_segment(segment_list(daf), 3, eid), time)
                end
            else
                if lid == 4
                    return $(fun1)(daf, get_segment(segment_list(daf), 4, eid), time)
                elseif lid == 5
                    return $(fun1)(daf, get_segment(segment_list(daf), 5, eid), time)
                else 
                    return $(fun1)(daf, get_segment(segment_list(daf), 6, eid), time)
                end
            end
        end

        # TODO: write custom transform function for PCKs which have a limited subset of types

    end

end