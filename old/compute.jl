@inline function findsegment(segments, origin::Int, target::Int)
    if !(origin in keys(segments) || target in keys(segments))
        throw(ArgumentError("No segment '$origin'->'$target' available."))
    end
    sign = 1.0
    if target < origin
        origin, target = target, origin
        sign = -1.0
    end
    segments[origin][target], sign
end

function unsafe_position!(pos, spk::SPK, tdb::Float64, tdb2::Float64, from::Int, to::Int)
    seg, sign = findsegment(spk.segments, from, to)
    unsafe_position!(pos, spk, seg, sign, tdb, tdb2)
end

function unsafe_velocity!(vel, spk::SPK, tdb::Float64, tdb2::Float64, from::Int, to::Int)
    seg, sign = findsegment(spk.segments, from, to)
    unsafe_velocity!(vel, spk, seg, sign, tdb, tdb2)
end

function unsafe_state!(pos, vel, spk::SPK, tdb::Float64, tdb2::Float64, from::Int, to::Int)
    seg, sign = findsegment(spk.segments, from, to)
    unsafe_state!(pos, vel, spk, seg, sign, tdb, tdb2)
end

function position!(pos, eph::SPK, tdb::Float64, tdb2::Float64, from::Int, to::Int)
    pos .= 0.0 
    path = get_nodes(eph.graph, from, to)
    for (origin, target) in zip(path[1:end-1], path[2:end])
        unsafe_position!(pos, eph::SPK, tdb::Float64, tdb2::Float64, origin::Int, target::Int)
    end
    pos 
end

function velocity!(vel, eph::SPK, tdb::Float64, tdb2::Float64, from::Int, to::Int)
    vel .= 0.0 
    path = get_nodes(spk.graph, from, to)
    for (origin, target) in zip(path[1:end-1], path[2:end])
        unsafe_velocity!(vel, eph::SPK, tdb::Float64, tdb2::Float64, origin::Int, target::Int)
    end
    vel 
end

function state!(pos, vel, eph::SPK, tdb::Float64, tdb2::Float64, from::Int, to::Int)
    pos .= 0.0
    vel .= 0.0 
    path = get_nodes(eph.graph, from, to)
    for (origin, target) in zip(path[1:end-1], path[2:end])
        unsafe_state!(pos, vel, eph::SPK, tdb::Float64, tdb2::Float64, origin::Int, target::Int)
    end
    pos, vel 
end
