
""" 
    SPKSegmentHeader20(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 20
"""
function SPKSegmentHeader20(daf::DAF, desc::DAFSegmentDescriptor)

    i0 = 8*(final_address(desc)-7)

    # Length and time scales in km and seconds
    dscale = get_float(array(daf), i0, endian(daf))
    tscale = get_float(array(daf), i0 + 8, endian(daf))

    # Everything is now transformed in seconds past J2000 and kms 
    DJ2000, D2S = 2451545, 86400

    # Integer and fractional part of the initial juliad date
    initjd = get_float(array(daf), i0 + 16, endian(daf))
    initfr = get_float(array(daf), i0 + 24, endian(daf))

    # Start epoch and interval length (in seconds)
    tstart = ((initjd - DJ2000) + initfr)*D2S
    tlen = get_float(array(daf), i0 + 32, endian(daf))*D2S

    # Number of elements in each record
    rsize = Int(get_float(array(daf), i0 + 40, endian(daf)))

    # Byte size of each logical record
    recsize = 8*rsize 

    # Polynomial degree 
    order = (rsize - 3) ÷ 3 - 1

    # Polynomial group size (number of coefficients required for the interpolation)
    N = order + 1

    # Number of records 
    n = Int(get_float(array(daf), i0 + 48, endian(daf)))

    # Initial segment address 
    iaa = initial_address(desc)

    SPKSegmentHeader20(dscale, tscale, tstart, tlen, recsize, order, N, n, iaa)
end

""" 
    SPKSegmentCache20(spkhead::SPKSegmentHeader2)

Initialise the cache for an SPK segment of type 20.
"""
function SPKSegmentCache20(head::SPKSegmentHeader20) 
    SPKSegmentCache20(
        MVector(-1), 
        MVector(0.0, 0.0, 0.0), 
        zeros(3, head.N)
        )
end

""" 
    SPKSegmentType20(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 2.
"""
function SPKSegmentType20(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader20(daf, desc)
    caches = [SPKSegmentCache20(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType20(header, caches)

end

@inline spk_field(::SPKSegmentType20) = SPK_SEGMENTLIST_MAPPING[20]

function spk_vector3(daf::DAF, seg::SPKSegmentType20, time::Number) 

    head = header(seg)
    data = cache(seg)

    # TODO: implement me
    get_coefficients!(daf, head, data, time)


end


function spk_vector6(daf::DAF, seg::SPKSegmentType20, time::Number)

    head = header(seg)
    data = cache(seg)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index = find_logical_record(head, time)
    get_coefficients!(daf, head, data, index)

    # TODO: implement me
    x, y, z = 0.0, 0.0, 0.0

    # Normalise the time argument between [-1, 1]
    t = chebyshev_time(head, time, index)

    # Compute the Chebyshev polynomials
    cheb = chebyshev!(data, t, head.N)
    vel = interpol(data, cheb, 0, 1, 0, head.N)*head.dscale/head.tscale

    return SA[x, y, z, vel[1], vel[2], vel[3]]

end

function spk_vector9(daf::DAF, seg::SPKSegmentType20, time::Number)

    head = header(seg)
    data = cache(seg)

    # TODO: implement me

end

function spk_vector12(daf::DAF, seg::SPKSegmentType20, time::Number)

    head = header(seg)
    data = cache(seg)

    # TODO: implement me

end


"""
    find_logical_record(head::SPKSegmentHeader2, time::Number)
"""
function find_logical_record(head::SPKSegmentHeader20, time::Number)
    
    # The index is returned in 0-based notation (i.e., the first record has index 0)
    index = floor(Int, (time - head.tstart)/head.tlen)

    if index == head.n 
        # This happens only when the time equals the final segment time
        index -= 1 
    end 

    return index
end

"""
    get_coefficients!(daf::DAF, head, cache, index::Int)
"""
function get_coefficients!(
    daf::DAF, head::SPKSegmentHeader20, cache::SPKSegmentCache20, index::Int
    )

    # Check whether the coefficients for this record are already loaded
    index == cache.id[1] && return nothing
    cache.id[1] = index 

    # Address of desired logical record 
    k = 8*(head.iaa-1) + head.recsize*index

    # For type 20 we do not have repeated midpoint and radius values 
    @inbounds for j = 1:3
        
        for i = 1:head.N 
            cache.A[j, i] = get_float(array(daf), k, endian(daf))
            k += 8
        end

        cache.p[j] = get_float(array(daf), k, endian(daf))
        k += 8

    end

    nothing 

end

function chebyshev_time(head::SPKSegmentHeader20, time::Number, index::Int)
    tbeg = head.tstart + head.tlen*index
    hlen = head.tlen/2 
    return (time - tbeg)/hlen - 1
end



"""
    chebyshev!(cache::SPKSegmentCache2, t::Number, N::Int)
"""
function chebyshev!(cache::SPKSegmentCache20, t::Number, N::Int)

    x1 = zeros(max(3, N))

    @inbounds begin 
        x1[1] = 1
        x1[2] = t 
        x1[3] = 2t*t - 1
        
        for j = 4:N 
            x1[j] = 2*t*x1[j-1] - x1[j-2]
        end
    end

    return x1
end 

"""
    interpol(cache::SPKSegmentCache2, cheb, order::Int, scale::Number, offset::Int, N::Int)
"""
function interpol(cache::SPKSegmentCache20, cheb::AbstractVector{T}, order::Int, 
            scale::Number, offset::Int, N::Int) where T

    ax = 1 + offset
    ay = 2 + offset
    az = 3 + offset

    x, y, z = T(0), T(0), T(0)
    @inbounds @simd for i in order+1:N
        x += cheb[i]*cache.A[ax, i]
        y += cheb[i]*cache.A[ay, i]
        z += cheb[i]*cache.A[az, i]
    end

    return scale*SA[x, y, z]

end