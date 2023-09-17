

jEph.load(::Type{EphemerisProvider}, file::AbstractString) = EphemerisProvider([file])

function jEph.load(::Type{EphemerisProvider}, files::Vector{<:AbstractString})
    return EphemerisProvider(files)
end

function Base.convert(::Type{jEph.EphemPointRecord}, r::EphemRecordSPK)
    DJ2000 = 2451545
    return jEph.EphemPointRecord(
        r.target, r.center, r.t_start[1]/86400 + DJ2000, r.t_end[end]/86400 + DJ2000, r.axes
    )
end

function Base.convert(::Type{jEph.EphemAxesRecord}, r::EphemRecordPCK)
    DJ2000 = 2451545
    return jEph.EphemAxesRecord(
        r.target, r.t_start[1]/86400 + DJ2000, r.t_end[end]/86400 + DJ2000, r.center
    )
end 

"""
    ephem_position_records(eph::EphemerisProvider)

Get an array of `EphemPointRecord`, providing detailed informations on the 
position content of the ephemeris file.
"""
function jEph.ephem_position_records(eph::EphemerisProvider)
    convert.(jEph.EphemPointRecord, ephem_spk_records(eph))
end

"""
    ephem_available_points(eph::EphemerisProvider)

Return a list of NAIFIds representing bodies with available ephemeris data. 
"""
jEph.ephem_available_points(eph::EphemerisProvider) = ephem_get_points(eph)

"""
    ephem_orient_records(eph::EphemerisProvider)

Get an array of `EphemAxesRecord`, providing detailed 
informations on the orientation content of the ephemeris file.
"""
function jEph.ephem_orient_records(eph::EphemerisProvider)
    convert.(jEph.EphemAxesRecord, ephem_pck_records(eph))
end

"""
    ephem_available_points(eph::EphemerisProvider)

Return a list of Frame IDs representing axes with available orientation data. 
"""
jEph.ephem_available_axes(eph::EphemerisProvider) = ephem_get_axes(eph)


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

    ft, lt, ct = analyse_timespan([ephem_spk_records(eph)..., ephem_pck_records(eph)...])
    
    DJ2000 = 2451544.5 # TODO: why do i have to put .5? 
    return ft/86400 + DJ2000, lt/86400 + DJ2000, ct

end

"""
    ephem_timescale(eph::EphemerisProvider) 

Retrieve a timescale ID associated with the ephemeris handler `eph`. 
It returns 1 for Barycentric Dynamical Time (TDB) and 2 for Barycentric Coordinate Time (TCB).

!!! warning 
    An error is thrown if the timescale is neither TDB nor TCB.
"""
jEph.ephem_timescale(eph::EphemerisProvider) = ephem_timescale_id(eph)

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

    # Transform time in seconds since J2000.0
    tsec = 86400*((jd0 + time) - 2451545)

    if order < 2 
        if order == 0
            res[1:3] .= ephem_vector3(eph, center, target, tsec)
        elseif order == 1
            res[1:6] .= ephem_vector6(eph, center, target, tsec)
        end
    else 
        if order == 2 
            res[1:9] .= ephem_vector9(eph, center, target, tsec)
        elseif order == 3 
            res[1:12] .= ephem_vector12(eph, center, target, tsec)
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


"""
    ephem_orient!(res, eph, jd0, time, target, center, order)

Interpolate the orientation and its derivatives up to `order` for the `target` body at the 
time `jd0` + `time`, expressed as a Julian Date. This function reads the ephemeris files 
associated to `eph` and stores the results to `res`.

The returned array `res` must be large enough to store the results. The size of this array 
must be equal to 3*order: 

- res[1:3] contain the angles 
- res[4:6] contain the 1st derivative  for order ≥ 1 
- res[7:9] contain the 2nd derivative for order ≥ 2
- res[10:12] contain the 3rd derivative for order ≥ 3

The values stores in `res` are always returned in rad, rad/s, rad/s², rad/s³

"""
function jEph.ephem_orient!(
    res, eph::EphemerisProvider, jd0::Number, time::Number, 
    target::Int, center::Int, order::Int
)

    # Transform time in seconds since J2000.0
    tsec = 86400*((jd0 + time) - 2451545)

    if order < 2 
        if order == 0
            res[1:3] .= ephem_rotation3(eph, center, target, tsec)
        elseif order == 1
            res[1:6] .= ephem_rotation6(eph, center, target, tsec)
        end
    else 
        if order == 2 
            res[1:9] .= ephem_rotation9(eph, center, target, tsec)
        elseif order == 3 
            res[1:12] .= ephem_rotation12(eph, center, target, tsec)
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