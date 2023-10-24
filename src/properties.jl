
export  ephem_timescale_id,
        ephem_spk_timespan,         
        ephem_pck_timespan, 
        ephem_spk_records, 
        ephem_pck_records,
        ephem_get_points,
        ephem_get_axes

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
    ephem_timescale_id(eph::EphemerisProvider)

Retrieve a timescale ID associated with the ephemeris handler `eph`. 
It returns 1 for Barycentric Dynamical Time (TDB) and 2 for Barycentric Coordinate Time (TCB).

!!! warning 
    Ephemeris providers with mixed timescales are not supported. An error is thrown if in 
    the ephemeris handler some segments are defined in TDB and some other segments in TCB.
"""
function ephem_timescale_id(eph::EphemerisProvider)

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
                        jEph.EphemerisError(
                            "Both TCB and TDB timescales are present in the same provider."
                        )
                    )
                end
            else
                if timescale_id == -1 || timescale_id == 2
                    timescale_id = 2
                else 
                    throw(
                        jEph.EphemerisError(
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
    ephem_get_points(eph::EphemerisProvider)

Return a list of NAIFIds representing bodies with available ephemeris data. 
"""
ephem_get_points(eph::EphemerisProvider) = Int.(keys(spk_links(eph)))

"""
    ephem_get_axes(eph::EphemerisProvider)

Return a list of Frame IDs representing axes with available orientation data. 
"""
ephem_get_axes(eph::EphemerisProvider) = Int.(keys(pck_links(eph)))

"""
    AbstractEphemRecord

Abstract type for ephemeris segment records.
"""
abstract type AbstractEphemRecord end 

"""
    initial_times(record::AbstractEphemRecord)

Recover the initial times of each sub-window in which the ephemeris data of `record` is 
defined, expressed in TDB seconds since J2000
"""
function initial_times(::T) where {T <: AbstractEphemRecord}
    throw(ErrorException("`initial_times` shall be implemented for type $T."))
end

"""
    final_times(record::AbstractEphemRecord)

Recover the final times of each sub-window in which the ephemeris data of `record` is 
defined, expressed in TDB seconds since J2000
"""
function final_times(::T) where {T <: AbstractEphemRecord}
    throw(ErrorException("`final_times` shall be implemented for type $T."))
end

"""
    EphemRecordSPK <: AbstractEphemRecord

Store the SPK metadata relative to a given (target, center) objects pair.

### Fields 
- `target` -- `Int` target object NAIF ID
- `center` -- `Int` center object NAIF ID
- `axes` -- `Int` reference axes ID 
- `t_start` -- start times of each sub-window, in TDB seconds since J2000
- `t_end` -- final times of each sub-window, in TDB seconds since J2000

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

Return a list of `EphemRecordSPK` storing metadata relative to each (target, center)
object pairs in the loaded SPK kernels. 
"""
function ephem_spk_records(eph::EphemerisProvider)
    
    # Retrieve the list of all the SPK segment descriptors
    spk_desc = DAFSegmentDescriptor[]
    for file in get_daf(eph)
        if is_spk(file)
            push!(spk_desc, descriptors(file)...)   
        end 
    end

    # SPK records vector
    records = EphemRecordSPK[]

    body_parsed = Tuple{Int, Int}[]
    axes_parsed = Int[]

    for desc in spk_desc 
        tid = target(desc)
        cid = center(desc)
        aid = axes(desc)

        # Retrieve all the descriptors associated to the same pair of points
        if !((tid, cid) in body_parsed || (cid, tid) in body_parsed)
            
            # Update the list of parsed points and the first axes they are defined on
            push!(body_parsed, (tid, cid))
            push!(axes_parsed, aid)

            desclist = DAFSegmentDescriptor[]
            for desc_k in spk_desc 
                if ((target(desc_k) == tid && center(desc_k) == cid) || 
                    (center(desc_k) == tid && target(desc_k) == cid))

                    push!(desclist, desc_k)
                
                end
            end

            if length(unique(axes.(desclist))) != 1
                @warn "The ephemeris data from $cid to $tid is defined on different axes."*
                      "Only axes with ID $aid are shown in the records."
            end
    
            t_start, t_end = get_segment_boundaries(desclist)
            push!(records, EphemRecordSPK(tid, cid, axes(desclist[1]), t_start, t_end))

        end
    end

    return records
    
end

"""
    EphemRecordPCK <: AbstractEphemRecord

Store the PCK metadata relative to a given (target, center) axes pair.

### Fields 
- `target` -- `Int` target axes NAIF ID
- `center` -- `Int` center axes NAIF ID
- `t_start` -- start times of each sub-window, in TDB seconds since J2000
- `t_end` -- final times of each sub-window, in TDB seconds since J2000
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

Return a list of `EphemRecordPCK` storing metadata relative to each (target, center)
axes pairs in the loaded PCK kernels. 
"""
function ephem_pck_records(eph::EphemerisProvider)

    # Retrieve the list of all the SPK segment descriptors
    pck_desc = DAFSegmentDescriptor[]
    for file in get_daf(eph)
        if is_pck(file)
            push!(pck_desc, descriptors(file)...)   
        end 
    end

    # SPK records vector
    records = EphemRecordPCK[]
    axes_parsed = Tuple{Int, Int}[]

    for desc in pck_desc 
        tid = target(desc)
        cid = center(desc)

        # Retrieve all the descriptors associated to the same pair of points
        if !((tid, cid) in axes_parsed || (cid, tid) in axes_parsed)
            
            # Update the list of parsed points and the first axes they are defined on
            push!(axes_parsed, (tid, cid))

            desclist = DAFSegmentDescriptor[]
            for desc_k in pck_desc 
                if ((target(desc_k) == tid && center(desc_k) == cid) || 
                    (center(desc_k) == tid && target(desc_k) == cid))

                    push!(desclist, desc_k)
                
                end
            end
    
            t_start, t_end = get_segment_boundaries(desclist)
            push!(records, EphemRecordPCK(tid, cid, t_start, t_end))

        end
    end

    return records
end


"""
    get_segment_boundaries(desclist::Vector{DAFSegmentDescriptor})

Parse all the segment descriptors of a given (center, target) pair and return a set of 
initial and final times, in TDB seconds since J2000, representing all the time 
sub-windows in which the ephemeris data for this pair is defined. 
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
            c = findlast(x -> x <= te, t_start)

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
            c = findfirst(x -> x >= ts, t_end)

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

        # if a == b the segment is discarded because 
        # its time interval is already covered.

    end

    return t_start, t_end

end

"""
    ephem_spk_timespan(eph::EphemerisProvider)

Return the minimum and maximum time available in the SPK kernels loaded within `eph`, in 
TDB seconds since J2000, together with a continuity parameter defined as follows: 

- **0** no SPK data is available.
- **1** the quantities of all bodies are available for any time between the first and last time.
- **2** the quantities of some bodies are available on discontinuous time intervals between the 
    first and last time.
- **3** the quantities of each body are available on a continuous time interval between the first 
    and the last time, but not available for any time between the first and last time.


### References 
- [CALCEPH](https://www.imcce.fr/inpop/calceph) C++ library

### See Also 
See also [`ephem_pck_timespan`](@ref).

"""
ephem_spk_timespan(eph::EphemerisProvider) = analyse_timespan(ephem_spk_records(eph))

"""
    ephem_pck_timespan(eph::EphemerisProvider)

Return the minimum and maximum time available in the PCK kernels loaded within `eph`, in 
TDB seconds since J2000, together with a continuity parameter defined as follows: 

- **0** no PCK data is available.
- **1** the quantities of all axes are available for any time between the first and last time.
- **2** the quantities of some axes are available on discontinuous time intervals between the 
    first and last time.
- **3** the quantities of each axis are available on a continuous time interval between the first 
    and the last time, but not available for any time between the first and last time.

### References 
- [CALCEPH](https://www.imcce.fr/inpop/calceph) C++ library

### See Also 
See also [`ephem_spk_timespan`](@ref).
"""
ephem_pck_timespan(eph::EphemerisProvider) = analyse_timespan(ephem_pck_records(eph))

"""
    analyse_timespan(records)

Analyse a set of [`AbstractEphemRecord`](@ref), returning the minimum and maximum covered 
times, in TDB seconds since J2000, together with a continuity parameter.

### References 
- [CALCEPH](https://www.imcce.fr/inpop/calceph) C++ library

### See Also 
See also [`ephem_spk_timespan`](@ref) and [`ephem_pck_timespan`](@ref).
"""
function analyse_timespan(records::Vector{<:AbstractEphemRecord})
    firsttime, lasttime = Inf, -Inf
    continuity = 0

    for record in records
        
        t_start = initial_times(record)
        t_end = final_times(record)

        # Update time boundaries 
        firsttime = min(firsttime, t_start[1])
        lasttime = max(lasttime, t_end[end])

        continuity = length(t_start) > 1 ? 2 : max(continuity, 1)

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