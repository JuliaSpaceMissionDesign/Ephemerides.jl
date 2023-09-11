
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
    ncomp = desc.segtype == 2 ? 3 : 6

    # Polynomial groupsize (number of coefficients required for the interpolation)
    N = (rsize - 2) ÷ ncomp

    # Polynomial order
    order = N - 1

    # Number of records
    n = Int(get_float(array(daf), i0+24, endian(daf)))

    # Initial segemtn address 
    iaa = initial_address(desc)

    SPKSegmentHeader2(tstart, tlen, order, N, n, recsize, ncomp, iaa)

end

""" 
    SPKSegmentCache2(spkhead::SPKSegmentHeader2)

Initialise the cache for an SPK segment of type 2 and 3.
"""
function SPKSegmentCache2(spkhead::SPKSegmentHeader2) 
    SPKSegmentCache2(
        zeros(spkhead.ncomp, spkhead.N), 
        MVector(0.0, 0.0, 0.0),
        DiffCache(zeros(max(3, spkhead.N))), 
        DiffCache(zeros(max(3, spkhead.N))),
        MVector(-1)
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

function spk_vector3(daf::DAF, seg::Union{SPKSegmentType2, SPKSegmentType3}, time::Number) 

    head = header(seg)
    data = cache(seg)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, head, data, time)

    # Normalise the time argument between [-1, 1]
    t = chebyshev_time(data, time)
    
    # Compute the Chebyshev polynomials
    chebyshev!(data, t, head.N)

    # Interpolate the body position
    return interpol(data, get_tmp(data.x1, t), 0, 1, 0)

end


function spk_vector6(daf::DAF, seg::SPKSegmentType2, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, head, data, time)
    
    # Normalise the time argument between [-1, 1]
    t = chebyshev_time(data, time)
    
    # Compute the Chebyshev polynomials
    chebyshev!(data, t, head.N)
    pos = interpol(data, get_tmp(data.x1, t), 0, 1, 0)

    # Compute 1st derivatives of Chebyshev polynomials 
    @inbounds scale = data.p[3]
    ∂chebyshev!(data, t, head.N)
    vel = interpol(data, get_tmp(data.x2, t), 1, scale, 0)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3]
    ]

end

function spk_vector9(daf::DAF, seg::SPKSegmentType2, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, head, data, time)

    # Normalise the time argument between [-1, 1]
    t = chebyshev_time(data, time)

    # Compute the Chebyshev polynomials
    chebyshev!(data, t, head.N)
    pos = interpol(data, get_tmp(data.x1, t), 0, 1, 0)

    # Compute 1st derivatives of Chebyshev polynomials 
    @inbounds scale = data.p[3]
    ∂chebyshev!(data, t, head.N)
    vel = interpol(data, get_tmp(data.x2, t), 1, scale, 0)

    # Compute 2nd derivative of Chebyshev polynomials 
    ∂²chebyshev!(data, t, head.N)
    acc = interpol(data, get_tmp(data.x1, t), 2, scale*scale, 0)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3], 
        acc[1], acc[2], acc[3]
    ]

end

function spk_vector12(daf::DAF, seg::SPKSegmentType2, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, head, data, time)

    # Normalise the time argument between [-1, 1]
    t = chebyshev_time(data, time)
    
    # Compute the Chebyshev polynomials
    chebyshev!(data, t, head.N)
    pos = interpol(data, get_tmp(data.x1, t), 0, 1, 0)

    # Compute 1st derivatives of Chebyshev polynomials 
    @inbounds scale = data.p[3]
    ∂chebyshev!(data, t, head.N)
    vel = interpol(data, get_tmp(data.x2, t), 1, scale, 0)

    # Compute 2nd derivative of Chebyshev polynomials 
    ∂²chebyshev!(data, t, head.N)
    acc = interpol(data, get_tmp(data.x1, t), 2, scale*scale, 0)

    # Compute 3rd derivative of Chebyshev polynomials 
    ∂³chebyshev!(data, t, head.N)
    jer = interpol(data, get_tmp(data.x2, t), 3, scale*scale*scale, 0)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3], 
        acc[1], acc[2], acc[3],
        jer[1], jer[2], jer[3]
    ]

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
    get_coefficients!(daf::DAF, head, cache, index::Int)
"""
function get_coefficients!(daf::DAF, head::SPKSegmentHeader2, cache::SPKSegmentCache2, 
            time::Number)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index = find_logical_record(head, time)

    # Check whether the coefficients for this record are already loaded
    index == cache.id[1] && return nothing
    cache.id[1] = index 

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

@inline function chebyshev_time(cache::SPKSegmentCache2, time::Number)
    @inbounds return (time - cache.p[1])/cache.p[2]
end

"""
    interpol(cache::SPKSegmentCache2, cheb, order::Int, scale::Number, offset::Int)
"""
function interpol(cache::SPKSegmentCache2, cheb::AbstractVector{T}, order::Int, 
            scale::Number, offset::Int) where T

    len = length(cheb)

    ax = 1 + offset
    ay = 2 + offset
    az = 3 + offset

    x, y, z = T(0), T(0), T(0)
    @inbounds @simd for i in order+1:len
        x += cheb[i]*cache.A[ax, i]
        y += cheb[i]*cache.A[ay, i]
        z += cheb[i]*cache.A[az, i]
    end

    return scale*SA[x, y, z]

end

"""
    chebyshev!(cache::SPKSegmentCache2, t::Number, N::Int)
"""
function chebyshev!(cache::SPKSegmentCache2, t::Number, N::Int)

    x1 = get_tmp(cache.x1, t)

    @inbounds begin 
        x1[1] = 1
        x1[2] = t 
        x1[3] = 2t*t - 1
        
        for j = 4:N 
            x1[j] = 2*t*x1[j-1] - x1[j-2]
        end
    end

    nothing
end 

"""
    ∂chebyshev!(cache::SPKSegmentCache2, t::Number, N::Int)
"""
function ∂chebyshev!(cache::SPKSegmentCache2, t::Number, N::Int)

    x1 = get_tmp(cache.x1, t)
    x2 = get_tmp(cache.x2, t)

    @inbounds begin 

        x2[1] = 0
        x2[2] = 1
        x2[3] = 4*t 

        for i = 4:N
            x2[i] = 2*t*x2[i-1] + 2*x1[i-1] - x2[i-2]
        end
    end 

    nothing 
end 

"""
    ∂²chebyshev!(cache::SPKSegmentCache2, t::Number, N::Int)
"""
function ∂²chebyshev!(cache::SPKSegmentCache2, t::Number, N::Int)

    x1 = get_tmp(cache.x1, t)
    x2 = get_tmp(cache.x2, t)

    @inbounds begin 

        x1[1] = 0
        x1[2] = 0
        x1[3] = 4 

        for i = 4:N 
            x1[i] = 4*x2[i-1] + 2*t*x1[i-1] - x1[i-2]
        end

    end

    nothing
end

"""
    ∂³chebyshev!(cache::SPKSegmentCache2, t::Number, N::Int)
"""
function ∂³chebyshev!(cache::SPKSegmentCache2, t::Number, N::Int)

    x1 = get_tmp(cache.x1, t)
    x2 = get_tmp(cache.x2, t)

    @inbounds begin 

        x2[1] = 0
        x2[2] = 0
        x2[3] = 0 

        for i = 4:N 
            x2[i] = 6*x1[i-1] + 2*t*x2[i-1] - x2[i-2]
        end

    end

    nothing
end

