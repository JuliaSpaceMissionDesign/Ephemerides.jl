

jEph.load(::Type{EphemerisProvider}, file::AbstractString) = EphemerisProvider([file])

function jEph.load(::Type{EphemerisProvider}, files::Vector{<:AbstractString})
    return EphemerisProvider(files)
end


"""
    ephem_position_records(eph::EphemerisProvider)

Get an array of [`jEph.EphemPointRecord`](@ref), providing detailed informations on the 
position content of the ephemeris file.
"""
function jEph.ephem_position_records(eph::EphemerisProvider)
    # TODO: implement
    throw(ErrorException("NotImplementedError"))
end

"""
    ephem_available_points(eph::EphemerisProvider)

Return a list of NAIFIds representing bodies with available ephemeris data. 
"""
function jEph.ephem_available_points(eph::EphemerisProvider)
    rec = jEph.ephem_position_records(eph)
    tids = map(x -> x.target, rec)
    cids = map(x -> x.center, rec)

    return unique(Int64[tids..., cids...])
end

"""
    ephem_orient_records(eph::EphemerisProvider)

Get an array of [`jEph.EphemAxesRecord`](@ref), providing detailed 
informations on the orientation content of the ephemeris file.
"""
function jEph.ephem_orient_records(eph::EphemerisProvider)
    # TODO: implement
    throw(ErrorException("NotImplementedError"))
end

"""
    ephem_available_points(eph::EphemerisProvider)

Return a list of Frame IDs representing axes with available orientation data. 
"""
function jEph.ephem_available_axes(eph::EphemerisProvider)
    rec = jEph.ephem_orient_records(eph)

    tids = map(x -> x.target, rec)
    cids = map(x -> x.axes, rec)

    return unique(Int64[tids..., cids...])
end

"""
    ephem_timespan(eph::EphemerisProvider)

Returns the first and last time available in the ephemeris file associated to 
an ephemeris file. It returns a tuple containing:

- `firsttime` -- Julian date of the first time.
- `lasttime` -- Julian date of the last time.
- `continuous` -- Information about the availability of the quantities over the 
                  time span. It equals:
                
    - `1`: if the quantities of all bodies are available for any time between 
        the first and last time.
    - `2`: if the quantities of some bodies are available on discontinuous time 
        intervals between the first and last time.
    - `3`: if the quantities of each body are available on a continuous time 
        interval between the first and last time, but not available for any 
        time between the first and last time.

"""
function jEph.ephem_timespan(eph::EphemerisProvider)
    # TODO: implement
    throw(ErrorException("NotImplemetedError"))
end

"""
    ephem_timescale(eph::EphemerisProvider) 

Retrieve a timescale ID associated with the ephemeris handler `eph`. 
It returns 1 for Barycentric Dynamical Time (TDB) and 2 for Barycentric Coordinate Time (TCB).

!!! warning 
    An error is thrown if the timescale is neither TDB nor TCB.
"""
function jEph.ephem_timescale(eph::EphemerisProvider)
    # TODO: implement

    if tsid == 1 || tsid == 2
        return tsid
    else
        throw(jEph.EphemerisError("unknown time scale identifier: $tsid"))
    end
end

"""
    ephem_compute!(res, eph, jd0, time, target, center, order)

Interpolate the position and/or its derivatives up to `order` for one body `target` relative 
to another `center` at the time `jd0` + `time`, expressed as a Julian Date. This function reads 
the ephemeris files associated to `eph` and stores the results to `res`.

The returned array `res` must be large enough to store the results. The size of this array 
must be equal to 3*order: 

- res[1:3] contain the position (x, y, z) and is always valid 
- res[4:6] contain the velocity (dx/dt, dy/dt, dz/dt) for order ≥ 1 
- res[7:9] contain the acceleration (d²x/dt², d²y/dt², d²z/dt²) for order ≥ 2
- res[10:12] contain the jerk (d³x/dt³, d³y/dt³, d³z/dt³) for order ≥ 3

The values stores in `res` are always returned in km, km/s, km/s², km/s³

### See also 
See also [`ephem_orient!`](@ref)
"""
function jEph.ephem_compute!(
    res,
    eph::EphemerisProvider,
    jd0::Number,
    time::Number,
    target::Int,
    center::Int,
    order::Int,
)

    if order < 2 
        if order == 0
            res[1:3] .= ephem_vector3(eph, center, target, jd0+time)
        elseif order == 1
            res[1:6] .= ephem_vector6(eph, center, target, jd0+time)
        end
    else 
        if order == 2 
            res[1:9] .= ephem_vector9(eph, center, target, jd0+time)
        elseif order == 3 
            res[1:12] .= ephem_vector12(eph, center, target, jd0+time)
        else
            throw(
                jEph.EphemerisError(
                    "the maximum accepted order is 3."
                )
            )
        end
    end

    return nothing
end