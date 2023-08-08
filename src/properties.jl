
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

"""
    ephem_available_points(eph::EphemerisProvider)

Return a list of NAIFIds representing bodies with available ephemeris data. 
"""
ephem_available_points(eph::EphemerisProvider) = Int.(keys(spk_links(eph)))

"""
    ephem_available_axes(eph::EphemerisProvider)

Return a list of Frame IDs representing axes with available orientation data. 
"""
ephem_available_axes(eph::EphemerisProvider) = Int.(keys(pck_links(eph)))

"""
    AbstractEphemRecord
"""
abstract type AbstractEphemRecord end 

"""
    initial_times(record::AbstractEphemRecord)
"""
function initial_times(::T) where {T <: AbstractEphemRecord}
    throw(ErrorException("`initial_times` shall be implemented for type $T."))
end

"""
    final_times(record::AbstractEphemRecord)
"""
function final_times(::T) where {T <: AbstractEphemRecord}
    throw(ErrorException("`final_times` shall be implemented for type $T."))
end

"""
    EphemRecordSPK <: AbstractEphemRecord

Store the SPK metadata relative to a given (target, center) objects pair.

### Fields 
- `target` -- 
- `center` --
- `axes` --
- `t_start` --
- `t_end` --

"""
struct EphemRecordSPK <: AbstractEphemRecord
    target::Int 
    center::Int
    axes::Int 
    t_start::Vector{Float64}
    t_end::Vector{Float64}
end

@inline initial_times(record::EphemRecordSPK) = record.t_start 
@inline final_times(record::EphemRecordSPK) = record.t_end 


"""
    ephem_spk_records(eph::EphemerisProvider)

"""
function ephem_spk_records(eph::EphemerisProvider)
    
    records = EphemRecordSPK[]
    parsed = Tuple{Int, Int}[]

    links = spk_links(eph)
    for center in sort(ephem_available_points(eph), rev=true)
        for target in sort(Int.(keys(links[center])), rev=true)

            if !((center, target) in parsed || (target, center) in parsed)

                # Retrieve all the segment descriptors for this (target, center) pair
                desclist = descriptor.(links[center][target])

                if length(unique(axes.(desclist))) != 1 
                    throw(
                        EphemerisError(
                            "The ephemeris data from $center to $target is defined on "*
                            "different axes."
                        )
                    )
                end

                t_start, t_end = get_segment_boundaries(desclist)
            
                push!(records, EphemRecordSPK(target, center, axes(desclist[1]), t_start, t_end))
                push!(parsed, (center, target))
            end

        end
    end

    return records
    
end

"""
    EphemRecordPCK <: AbstractEphemRecord

Store the PCK metadata relative to a given (target, center) axes pair.

### Fields 
- `target` -- 
- `center` --
- `t_start` --
- `t_end` --
"""
struct EphemRecordPCK <: AbstractEphemRecord
    target::Int 
    center::Int
    t_start::Vector{Float64}
    t_end::Vector{Float64}
end

@inline initial_times(record::EphemRecordPCK) = record.t_start 
@inline final_times(record::EphemRecordPCK) = record.t_end 


"""
    ephem_pck_records(eph::EphemerisProvider)

"""
function ephem_pck_records(eph::EphemerisProvider)

    records = EphemRecordPCK[]
    parsed = Tuple{Int, Int}[]

    links = pck_links(eph)
    for center in sort(ephem_available_axes(eph), rev=true)
        for target in sort(Int.(keys(links[center])), rev=true)

            if !((center, target) in parsed || (target, center) in parsed)

                # Retrieve all the segment descriptors for this (target, center) pair
                desclist = descriptor.(links[center][target])
                t_start, t_end = get_segment_boundaries(desclist)
            
                push!(records, EphemRecordPCK(target, center, t_start, t_end))
                push!(parsed, (center, target))
            end

        end
    end

    return records
end


"""
    get_segment_boundaries(desclist::Vector{DAFSegmentDescriptor})
"""
function get_segment_boundaries(desclist::Vector{DAFSegmentDescriptor})

    t_start = Float64[]
    t_end = Float64[]

    for desc in desclist

        ts = initial_time(desc)
        te = final_time(desc)

        # TODO: what happens when the time boundaries match?
        a = findlast(x -> x < ts, t_start)
        b = findfirst(x -> x > te, t_end)

        if isnothing(a) && isnothing(b)
            t_start = [ts]
            t_end = [te]

        elseif isnothing(a) 
            c = findlast(x -> x < te, t_start)

            if isnothing(c)
                # The segment is before all the others
                insert!(t_start, 1, ts)
                insert!(t_end, 1, te)
            else 
                # Update t_start vector
                deleteat!(t_start, 1:c)
                insert!(t_start, 1, ts)
            
                # Update t_end vector
                deleteat!(t_end, 1:c-1)
                if c != b
                    t_end[1] = te
                end
            end

        elseif isnothing(b)
            c = findfirst(x -> x > ts, t_end)

            if isnothing(c)
                # The segment is after all the others
                push!(t_start, ts)
                push!(t_end, te)
            else
                # Update t_end vector 
                deleteat!(t_end, c:length(t_end))
                push!(t_end, te)

                # Update t_start vector 
                deleteat!(t_start, a+1:length(t_start))
                if c != a 
                    push!(t_start, ts)
                end

            end

        elseif a != b
            # Update t_start vector 
            deleteat!(t_start, a+1:c)
            d > a && insertat!(t_start, a+1, ts)
        
            # Update t_end vector 
            deleteat!(t_end, d:b-1)
            c < b && insertat!(t_end, c, te)

        end

    end

    return t_start, t_end

end


"""
    ephem_spk_timespan(eph::EphemerisProvider)
"""
ephem_spk_timespan(eph::EphemerisProvider) = analyse_timespan(ephem_spk_records(eph))


"""
    ephem_pck_timespan(eph::EphemerisProvider)
"""
ephem_pck_timespan(eph::EphemerisProvider) = analyse_timespan(ephem_pck_records(eph))


"""
    analyse_timespan(records)
"""
function analyse_timespan(records)
    firsttime, lasttime = Inf, -Inf
    continuity = 0

    for record in records
        
        t_start = initial_times(record)
        t_end = final_times(record)

        # Update time boundaries 
        firsttime = min(firsttime, t_start[1])
        lasttime = max(lasttime, t_end[end])

        continuity = length(t_start) > 1 ? 2 : continuity

    end

    # Check if all the records have the same timespan 
    if continuity == 1
        for record in records
            if (initial_times(record)[1] != firsttime || final_times(record)[end] != lasttime)
                continuity = 3 
            end
        end
    end

    return firsttime, lasttime, continuity

end