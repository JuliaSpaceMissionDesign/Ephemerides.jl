

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

function jEph.ephem_position_records(eph::EphemerisProvider)
    convert.(jEph.EphemPointRecord, ephem_spk_records(eph))
end

function jEph.ephem_orient_records(eph::EphemerisProvider)
    convert.(jEph.EphemAxesRecord, ephem_pck_records(eph))
end

jEph.ephem_available_points(eph::EphemerisProvider) = ephem_get_points(eph)
jEph.ephem_available_axes(eph::EphemerisProvider) = ephem_get_axes(eph)

function jEph.ephem_timespan(eph::EphemerisProvider)

    ft, lt, ct = analyse_timespan([ephem_spk_records(eph)..., ephem_pck_records(eph)...])
    
    DJ2000 = 2451545
    return ft/86400 + DJ2000, lt/86400 + DJ2000, ct

end

jEph.ephem_timescale(eph::EphemerisProvider) = ephem_timescale_id(eph)

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