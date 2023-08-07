

function ephem_vector3(eph::EphemerisProvider, from::Int, to::Int, time::Number)

    # TODO: check that the path actually exists
    for link in spk_links(eph)[from][to] 
        if initial_time(link) <= time <= final_time(link)   
            return factor(link)*ephem_vector3(get_daf(eph, file_id(link)), link, time)
        end
    end

    throw(
        jEph.EphemerisError(
            "ephemeris data for point with NAIFId $target with respect to point " *
            "$center is not available at JD $(jd0+time)"
        ),
    )

end

function ephem_vector3(daf::DAF, link::SPKLink, time::Number)

    # Retrieve list and element link IDs
    lid = list_id(link)
    eid = element_id(link)
    
    if lid == 1 
        return ephem_vector3(
            daf, get_segment(get_segment_list(daf), 1, eid), descriptor(link), time
        )
    else 
        return ephem_vector3(
            daf, get_segment(get_segment_list(daf), 2, eid), descriptor(link), time
        )
    end
        
end