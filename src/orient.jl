

function ephem_vector3(eph::EphemerisProvider, from::Int, to::Int, time::Number)

    # TODO: check that the path actually exists
    for link in eph.linktable[from][to] 
        if link.desc.tstart <= time <= link.desc.tend   
            return link.fct*ephem_vector3(eph.files[link.fid], link, time)
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

    if link.fid == 1 
        return ephem_vector3(daf, getfield(daf.seglist, 1)[link.lid], link.desc, time)
    else 
        return ephem_vector3(daf, getfield(daf.seglist, 2)[link.lid], link.desc, time)
    end
        
end