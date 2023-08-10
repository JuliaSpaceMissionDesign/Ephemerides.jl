
""" 
    SPKSegmentHeader2(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 2 and 3.
"""
function SPKSegmentHeader2(daf::DAF, desc::DAFSegmentDescriptor)

    i0 = 8*(final_address(desc)-4)

    tstart = get_float(array(daf), i0, endian(daf))
    tlen = get_float(array(daf), i0+8, endian(daf))

    # polynomial order 
    rsize = Int(get_float(array(daf), i0+16, endian(daf)))

    # Number of components 
    ncomp = desc.segtype == 2 ? 3 : 6

    # The order of the polynomial is actually = (order-1)
    order = (rsize - 2) ÷ ncomp

    # number of records
    n = Int(get_float(array(daf), i0+24, endian(daf)))

    # Byte size of each logical record
    recsize = 8*(order*ncomp + 2)

    # Chebyshev scale factor 
    scale = 2/tlen;

    # Initial segemtn address 
    iaa = initial_address(desc)

    SPKSegmentHeader2(tstart, tlen, order, n, recsize, ncomp, scale, iaa)

end

""" 
    SPKSegmentCache2(spkhead::SPKSegmentHeader2)

Initialise the cache for an SPK segment of type 2 and 3.
"""
function SPKSegmentCache2(spkhead::SPKSegmentHeader2) 
    SPKSegmentCache2(
        zeros(spkhead.ncomp, spkhead.order), 
        DiffCache(zeros(spkhead.order)), 
        DiffCache(zeros(spkhead.order)),
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


function spk_vector3(daf::DAF, seg::SPKSegmentType2, time::Number) 

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(header(seg), time)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, header(seg), cache(seg), index)
    
    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)

    # Interpolate the body position
    return interpol(cache(seg), get_tmp(cache(seg).x1, time), 0, 1, 0)
end


function spk_vector6(daf::DAF, seg::SPKSegmentType2, time::Number)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(header(seg), time)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, header(seg), cache(seg), index)

    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)
    pos = interpol(cache(seg), get_tmp(cache(seg).x1, time), 0, 1, 0)

    # Compute 1st derivatives of Chebyshev polynomials 
    ∂chebyshev!(cache(seg), t, seg.head.order)
    vel = interpol(cache(seg), get_tmp(cache(seg).x2, time), 1, seg.head.scale, 0)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3]
    ]

end

function spk_vector9(daf::DAF, seg::SPKSegmentType2, time::Number)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(header(seg), time)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, header(seg), cache(seg), index)

    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)
    pos = interpol(cache(seg), get_tmp(cache(seg).x1, time), 0, 1, 0)

    # Compute 1st derivatives of Chebyshev polynomials 
    scale = seg.head.scale 
    ∂chebyshev!(cache(seg), t, seg.head.order)
    vel = interpol(cache(seg), get_tmp(cache(seg).x2, time), 1, scale, 0)

    # Compute 2nd derivative of Chebyshev polynomials 
    scale *= scale
    ∂²chebyshev!(cache(seg), t, seg.head.order)
    acc = interpol(cache(seg), get_tmp(cache(seg).x1, time), 2, scale, 0)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3], 
        acc[1], acc[2], acc[3]
    ]

end

function spk_vector12(daf::DAF, seg::SPKSegmentType2, time::Number)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(header(seg), time)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, header(seg), cache(seg), index)
    
    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)
    pos = interpol(cache(seg), get_tmp(cache(seg).x1, time), 0, 1, 0)

    # Compute 1st derivatives of Chebyshev polynomials
    scale = seg.head.scale 
    ∂chebyshev!(cache(seg), t, seg.head.order)
    vel = interpol(cache(seg), get_tmp(cache(seg).x2, time), 1, scale, 0)

    # Compute 2nd derivative of Chebyshev polynomials 
    scale *= scale
    ∂²chebyshev!(cache(seg), t, seg.head.order)
    acc = interpol(cache(seg), get_tmp(cache(seg).x1, time), 2, scale, 0)

    # Compute 3rd derivative of Chebyshev polynomials 
    scale *= scale
    ∂³chebyshev!(cache(seg), t, seg.head.order)
    jer = interpol(cache(seg), get_tmp(cache(seg).x2, time), 3, scale, 0)

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

    idx, Δt = divrem(time - head.tstart, head.tlen) 
    index = round(Int, idx)

    if index == head.n 
        # This should only happen when time equals the final segment time
        index -= 1 
        Δt = head.tlen 
    end 

    # Time argument for the Chebyshev polynomials
    t = Δt*head.scale - 1 

    return index, t
end

"""
    get_coefficients!(daf::DAF, head, cache, index::Int)
"""
function get_coefficients!(daf::DAF, head::SPKSegmentHeader2, cache::SPKSegmentCache2, 
            index::Int)
        
    # Check whether the coefficients for this record are already loaded
    index == cache.id[1] && return nothing
    cache.id[1] = index 

    # Address of desired logical record 
    # (skipping mid and radius because they are all equal for SPK type 2)
    k = 8*(head.iaa-1) + head.recsize*index + 16

    # TODO: can we speed-up this part by casting the byte content into the array at once?
    @inbounds for j = 1:head.order 
        for i = 1:head.ncomp
            cache.A[i, j] = get_float(
                array(daf), k + 8*(j-1) + 8*(i-1)*head.order, endian(daf)
            )
        end
    end

    nothing

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
    chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)
"""
function chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)

    x1 = get_tmp(cache.x1, t)

    @inbounds begin 
        x1[1] = 1
        x1[2] = t 
        x1[3] = 2t*t - 1
        
        for i = 4:order 
            x1[i] = 2*t*x1[i-1] - x1[i-2]
        end
    end

    nothing
end 

"""
    ∂chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)
"""
function ∂chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)

    x1 = get_tmp(cache.x1, t)
    x2 = get_tmp(cache.x2, t)

    @inbounds begin 

        x2[1] = 0
        x2[2] = 1
        x2[3] = 4*t 

        for i = 4:order 
            x2[i] = 2*t*x2[i-1] + 2*x1[i-1] - x2[i-2]
        end
    end 

    nothing 
end 

"""
    ∂²chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)
"""
function ∂²chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)

    x1 = get_tmp(cache.x1, t)
    x2 = get_tmp(cache.x2, t)

    @inbounds begin 

        x1[1] = 0
        x1[2] = 0
        x1[3] = 4 

        for i = 4:order 
            x1[i] = 4*x2[i-1] + 2*t*x1[i-1] - x2[i-2]
        end

    end

    nothing
end

"""
    ∂³chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)
"""
function ∂³chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)

    x1 = get_tmp(cache.x1, t)
    x2 = get_tmp(cache.x2, t)

    @inbounds begin 

        x2[1] = 0
        x2[2] = 0
        x2[3] = 0 

        for i = 4:order 
            x2[i] = 6*x1[i-1] + 2*t*x2[i-1] - x2[i-2]
        end

    end

    nothing
end

