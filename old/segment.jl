function get_descriptor(record::Vector{UInt8}, little::Bool)
    firstsec, lastsec = reinterpret_getindex(Float64, record, (1, 9), little)
    target, center, frame, spktype, firstaddr, lastaddr =
            reinterpret_getindex(Int32, record, (17, 21, 25, 29, 33, 37), little)
    firstsec, lastsec, target, center, frame, spktype, firstaddr, lastaddr
end

# ----
# Type 2

function Type02Segment(daf::DAF, name::String, record::Vector{UInt8})
    firstsec, lastsec, target, center, frame, spktype, firstaddr, lastaddr = 
        get_descriptor(record, daf.little)
    
    i0 = lastaddr * SIZE_FLOAT64 - 4 * SIZE_FLOAT64 + 1
    init, intlen, rsize, n_records =
        reinterpret_getindex(Float64, daf.array, (i0, i0 + 8, i0 + 16, i0 + 24), daf.little)
    n_records = round(Int32, n_records)
    order = Int((rsize - 2) ÷ 3)
    Type02Segment(
        SPKDescriptor(
            name,
            firstsec,
            lastsec,
            jd(firstsec),
            jd(lastsec),
            target,
            center,
            frame,
            spktype,
            firstaddr,
            lastaddr
        ),
        firstaddr*SIZE_FLOAT64 - SIZE_FLOAT64 + 1,
        lastaddr*SIZE_FLOAT64 - SIZE_FLOAT64*4,
        init,
        intlen,
        round(Int32, rsize),
        n_records,
        order,
        -1,
        zeros(3, order),
        zeros(order),
        zeros(order),
    )
end

function parse_segment(daf::DAF, name::String, record::Vector{UInt8}, type::Val{2})
    Type02Segment(daf, name, record)
end

@inline function getrecordnum(seg::Type02Segment, tdb::Float64, tdb2::Float64)
    checkdate(seg, tdb+tdb2)
    secs = (seconds(tdb) - seg.initialsecond) + tdb2 * SECONDS_PER_DAY
    recordnum, frac = divrem(secs, seg.intlen)
    recordnum = round(Int, recordnum)
    if recordnum == seg.n_records
        recordnum -= 1
        frac = seg.intlen
    end
    recordnum, frac
end

@inline function update_cache!(spk::SPK, seg::Type02Segment, recordnum)
    components = 3
    seg.cached_record = recordnum
    # Drop the MID and RADIUS values
    first = seg.firstword + SIZE_FLOAT64 * seg.rsize * recordnum + SIZE_FLOAT64 * 2
    ptr = Ptr{Float64}(pointer(spk.daf.array, first))

    cache = unsafe_wrap(Array, Ptr{Float64}(ptr), (seg.order, components), own=false)
    if !spk.daf.little
        transpose!(seg.cache, ntoh.(copy(cache)))
    else
        transpose!(seg.cache, cache)
    end
end

@inline function chebyshev!(seg, frac)
    seg.x[1] = 1.0
    seg.x[2] = 2.0 * frac / seg.intlen - 1.0
    @inbounds for i = 3:seg.order
        seg.x[i] = 2.0 * seg.x[2] * seg.x[i-1] - seg.x[i-2]
    end
end

@inline function chebyshev_deriv!(seg)
    seg.t[2] = 1.0
    if seg.order > 2
        seg.t[3] = 4.0 * seg.x[2]
        @inbounds for i = 4:seg.order
            seg.t[i] = 2.0 * seg.x[2] * seg.t[i-1] - seg.t[i-2] +
                seg.x[i-1] + seg.x[i-1]
        end
    end
    seg.t .*= 2.0
    seg.t ./= seg.intlen
end

@inline function unsafe_position!(r, seg::Type02Segment, sign::Float64)
    @inbounds @simd for i = 1:3
        for j = 1:seg.order
            r[i] += sign * seg.cache[i, j] * seg.x[j]
        end
    end
    r
end

@inline function unsafe_position!(r, spk::SPK, seg::Type02Segment, sign::Float64, tdb::Float64, tdb2::Float64=0.0)
    recordnum, frac = getrecordnum(seg, tdb, tdb2)
    if recordnum != seg.cached_record
        update_cache!(spk, seg, recordnum)
    end
    chebyshev!(seg, frac)
    unsafe_position!(r, seg, sign)
end

@inline function unsafe_velocity!(v, seg::Type02Segment, sign::Float64)
    chebyshev_deriv!(seg)
    @inbounds @simd for i = 1:3
        for j = 1:seg.order
            v[i] += sign * seg.cache[i, j] * seg.t[j]
        end
    end
    v
end

@inline function unsafe_velocity!(v, spk::SPK, seg::Type02Segment, sign::Float64, tdb::Float64, tdb2::Float64=0.0)
    recordnum, frac = getrecordnum(seg, tdb, tdb2)
    if recordnum != seg.cached_record
        update_cache!(spk, seg, recordnum)
    end
    chebyshev!(seg, frac)
    unsafe_velocity!(v, seg, sign)
end

@inline function unsafe_state!(pos, vel, spk::SPK, seg::Type02Segment, sign::Float64, tdb::Float64, tdb2::Float64=0.0)
    recordnum, frac = getrecordnum(seg, tdb, tdb2)
    if recordnum != seg.cached_record
        update_cache!(spk, seg, recordnum)
    end
    chebyshev!(seg, frac)
    unsafe_position!(pos, seg, sign)
    unsafe_velocity!(vel, seg, sign)
    pos, vel
end
