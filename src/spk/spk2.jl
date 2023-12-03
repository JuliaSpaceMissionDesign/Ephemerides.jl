
""" 
    SPKSegmentHeader2(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 2 and 3.
"""
function SPKSegmentHeader2(daf::DAF, desc::DAFSegmentDescriptor)

    i0 = 8*(final_address(desc)-4)

    tstart = get_float(array(daf), i0, endian(daf))
    tlen = get_float(array(daf), i0+8, endian(daf))

    # Number of elements in each record
    rsize = Int(get_float(array(daf), i0+16, endian(daf)))

    # Byte size of each logical record
    recsize = 8*rsize

    # Number of components 
    ncomp = segment_type(desc) == 2 ? 3 : 6

    # Polynomial groupsize (number of coefficients required for the interpolation)
    N = (rsize - 2) ÷ ncomp

    # Polynomial order
    order = N - 1

    # Number of records
    n = Int(get_float(array(daf), i0+24, endian(daf)))

    # Initial segemtn address 
    iaa = initial_address(desc)

    SPKSegmentHeader2(tstart, tlen, order, N, n, recsize, ncomp, iaa, segment_type(desc))

end

""" 
    SPKSegmentCache2(head::SPKSegmentHeader2)

Initialise the cache for an SPK segment of type 2 and 3.
"""
function SPKSegmentCache2(head::SPKSegmentHeader2) 
    SPKSegmentCache2(
        zeros(head.ncomp, max(3, head.N)), 
        MVector(0.0, 0.0, 0.0),
        InterpCache{Float64}(4, max(3, head.N)),
        -1
    )
end

""" 
    SPKSegmentType2(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 2.
"""
function SPKSegmentType2(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader2(daf, desc)
    caches = [SPKSegmentCache2(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType2(header, caches)

end

@inline spk_field(::SPKSegmentType2) = SPK_SEGMENTLIST_MAPPING[2]

function spk_vector3(daf::DAF, seg::SPKSegmentType2, time::Number) 

    head = header(seg)
    data = cache(seg)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, head, data, time)

    # Normalise the time argument between [-1, 1]
    t = normalise_time(data, time)

    x, y, z = chebyshev(data.buff, data.A, t, 0, head.N)
    return SVector{3}(x, y, z)

end


function spk_vector6(daf::DAF, seg::SPKSegmentType2, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, head, data, time)
    
    # Normalise the time argument between [-1, 1]
    t = normalise_time(data, time)

    @inbounds if head.type == 2
        x, y, z, vx, vy, vz = ∂chebyshev(data.buff, data.A, t, 0, head.N, data.p[3])
    else 
        x, y, z = chebyshev(data.buff, data.A, t, 0, head.N)
        vx, vy, vz = chebyshev(data.buff, data.A, t, 3, head.N)
    end

    return SVector{6}(x, y, z, vx, vy, vz)

end

function spk_vector9(daf::DAF, seg::SPKSegmentType2, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, head, data, time)

    # Normalise the time argument between [-1, 1]
    t = normalise_time(data, time)

    @inbounds if head.type == 2
        x, y, z, vx, vy, vz, ax, ay, az = ∂²chebyshev(
            data.buff, data.A, t, 0, head.N, data.p[3]
        )

    else 
        x, y, z = chebyshev(data.buff, data.A, t, 0, head.N)
        vx, vy, vz, ax, ay, az = ∂chebyshev(data.buff, data.A, t, 3, head.N, data.p[3])
    end

    return SVector{9}(x, y, z, vx, vy, vz, ax, ay, az)

end

function spk_vector12(daf::DAF, seg::SPKSegmentType2, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, head, data, time)

    # Normalise the time argument between [-1, 1]
    t = normalise_time(data, time)

    @inbounds if head.type == 2
        x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz = ∂³chebyshev(
            data.buff, data.A, t, 0, head.N, data.p[3]
        )

    else 
        x, y, z = chebyshev(data.buff, data.A, t, 0, head.N)
        vx, vy, vz, ax, ay, az, jx, jy, jz = ∂²chebyshev(
            data.buff, data.A, t, 3, head.N, data.p[3]
        )
    end

    return SVector{12}(x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz)

end


"""
    find_logical_record(head::SPKSegmentHeader2, time::Number)
"""
function find_logical_record(head::SPKSegmentHeader2, time::Number)
    
    # The index is returned in 0-based notation (i.e., the first record has index 0)
    index = floor(Int, (time - head.tstart)/head.tlen)

    if index == head.n 
        # This happens only when the time equals the final segment time
        index -= 1 
    end 

    return index
end

"""
    get_coefficients!(daf::DAF, head::SPKSegmentHeader2, cache::SPKSegmentCache2, index::Int)
"""
function get_coefficients!(daf::DAF, head::SPKSegmentHeader2, cache::SPKSegmentCache2, 
            time::Number)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index = find_logical_record(head, time)

    # Check whether the coefficients for this record are already loaded
    index == cache.id && return nothing
    cache.id = index 

    # Address of desired logical record 
    k = 8*(head.iaa-1) + head.recsize*index

    # Retrieve record mid point and radius
    cache.p[1] = get_float(array(daf), k, endian(daf))
    cache.p[2] = get_float(array(daf), k + 8, endian(daf))
    cache.p[3] = 1/cache.p[2]

    k += 16

    # TODO: can we speed-up this part by casting the byte content into the array at once?
    @inbounds for j = 1:head.N 
        for i = 1:head.ncomp
            cache.A[i, j] = get_float(
                array(daf), k + 8*(j-1) + 8*(i-1)*head.N, endian(daf)
            )
        end
    end

    nothing

end

"""
    normalise_time(cache::SPKSegmentCache2, time::Number)

Transform `time` in an interval between [-1, 1] for compliance with Chebyshev polynomials.
"""
@inline function normalise_time(cache::SPKSegmentCache2, time::Number)
    @inbounds return (time - cache.p[1])/cache.p[2]
end
