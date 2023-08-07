
"""
    TDB_SEGMENTS

List of the SPK/PCK segment types for which the time argument is expressed in the TDB scale. 
"""
const TDB_SEGMENTS = (1, 2, 3, 5, 8, 9, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21)

"""
    TCB_SEGMENTS

List of the SPK/PCK segment types for which the time argument is expressed in the TCB scale.
"""
const TCB_SEGMENTS = (102, 103, 120)


"""
    eph_timescale(eph::EphemerisProvider)

Retrieve a timescale ID associated with the ephemeris handler `eph`. 
It returns 1 for Barycentric Dynamical Time (TDB) and 2 for Barycentric Coordinate Time (TCB).

!!! warning 
    Ephemeris providers with mixed timescales are not supported. An error is thrown if in 
    the ephemeris handler some segments are defined in TDB and some other segments in TCB.
"""
function ephem_timescale(eph::EphemerisProvider)

    timescale_id = -1 

    # Parse all the segment descriptors in all the DAFs contained 
    # in the ephemeris provider.
    for daf in get_daf(eph)
        for desc in descriptors(daf)

            if segment_type(desc) in TDB_SEGMENTS
                if timescale_id == -1 || timescale_id == 1
                    timescale_id = 1
                else 
                    throw(
                        EphemerisError(
                            "Both TCB and TDB timescales are present in the same provider."
                        )
                    )
                end
            else
                if timescale_id == -1 || timescale_id == 2
                    timescale_id = 2
                else 
                    throw(
                        EphemerisError(
                            "Both TCB and TDB timescales are present in the same provider."
                        )
                    )
                end
            end

        end
    end

    return timescale_id
    
end


function ephem_timespan(eph::EphemerisProvider)
    # TODO: 
end

function ephem_available_points(eph::EphemerisProvider)
    # TODO:
end

function ephem_available_axes(eph::EphemerisProvider)
    # TODO: 
end

function ephem_spk_records(eph::EphemerisProvider)
    # TODO:
end

function ephem_pck_records(eph::EphemerisProvider)
    # TODO: 
end