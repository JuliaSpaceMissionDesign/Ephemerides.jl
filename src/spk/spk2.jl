
""" 
    SPKSegmentHeader2(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 2.
"""
function SPKSegmentHeader2(daf::DAF, desc::DAFSegmentDescriptor)

    i0 = 8*(desc.faa-4)

    tstart = get_float(array(daf), i0, endian(daf))
    tlen = get_float(array(daf), i0+8, endian(daf))

    # polynomial order 
    rsize = Int(get_float(array(daf), i0+16, endian(daf)))

    # The order of the polynomial is actually = (order-1)
    order = (rsize - 2) ÷ 3

    # number of records
    n = Int(get_float(array(daf), i0+24, endian(daf)))

    # Number of components 
    ncomp = 3

    # Byte size of each logical record
    recsize = 8*(order*ncomp + 2)

    # Chebyshev scale factor 
    scale = 2/tlen;

    SPKSegmentHeader2(tstart, tlen, order, n, recsize, ncomp, scale)

end

""" 
    SPKSegmentCache2(spkhead::SPKSegmentHeader2)

Initialise the cache for an SPK segment of type 2.
"""
function SPKSegmentCache2(spkhead::SPKSegmentHeader2) 
    SPKSegmentCache2(
        zeros(3, spkhead.order), 
        zeros(spkhead.order), 
        zeros(spkhead.order), 
        MVector(-1)
    )
end

""" 
    SPKSegmentType2(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 2.
"""
function SPKSegmentType2(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    head = SPKSegmentHeader2(daf, desc)
    cache = SPKSegmentCache2(head)

    SPKSegmentType2(head, cache)

end

@inline spk_field(::SPKSegmentType2) = SPK_SEGMENTLIST_MAPPING[2]


function spk_vector3(daf::DAF, seg::SPKSegmentType2, desc::DAFSegmentDescriptor, time::Number) 

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(seg, time)
    
    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, seg, desc, index)

    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)

    # Interpolate the body position
    return interpol(seg.cache, cache(seg).x1, 0, 1)
end


function spk_vector6(daf::DAF, seg::SPKSegmentType2, desc::DAFSegmentDescriptor, time::Number)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(seg, time)
    
    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, seg, desc, index)

    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)
    pos = interpol(cache(seg), cache(seg).x1, 0, 1)

    # Compute 1st derivatives of Chebyshev polynomials 
    ∂chebyshev!(cache(seg), t, seg.head.order)
    vel = interpol(cache(seg), cache(seg).x2, 1, seg.head.scale)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3]
    ]

end

function spk_vector9(daf::DAF, seg::SPKSegmentType2, desc::DAFSegmentDescriptor, time::Number)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(seg, time)
    
    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, seg, desc, index)

    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)
    pos = interpol(cache(seg), cache(seg).x1, 0, 1)

    # Compute 1st derivatives of Chebyshev polynomials 
    scale = seg.head.scale 
    ∂chebyshev!(cache(seg), t, seg.head.order)
    vel = interpol(cache(seg), cache(seg).x2, 1, scale)

    # Compute 2nd derivative of Chebyshev polynomials 
    scale *= scale
    ∂²chebyshev!(cache(seg), t, seg.head.order)
    acc = interpol(cache(seg), cache(seg).x1, 2, scale)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3], 
        acc[1], acc[2], acc[3]
    ]

end

function spk_vector12(daf::DAF, seg::SPKSegmentType2, desc::DAFSegmentDescriptor, time::Number)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(seg, time)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, seg, desc, index)
    
    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)
    pos = interpol(cache(seg), cache(seg).x1, 0, 1)

    # Compute 1st derivatives of Chebyshev polynomials
    scale = seg.head.scale 
    ∂chebyshev!(cache(seg), t, seg.head.order)
    vel = interpol(cache(seg), cache(seg).x2, 1, scale)

    # Compute 2nd derivative of Chebyshev polynomials 
    scale *= scale
    ∂²chebyshev!(cache(seg), t, seg.head.order)
    acc = interpol(cache(seg), cache(seg).x1, 2, scale)

    # Compute 3rd derivative of Chebyshev polynomials 
    scale *= scale
    ∂³chebyshev!(cache(seg), t, seg.head.order)
    jer = interpol(cache(seg), cache(seg).x2, 3, scale)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3], 
        acc[1], acc[2], acc[3],
        jer[1], jer[2], jer[3]
    ]

end


"""
    find_logical_record(seg::SPKSegmentType2, time::Number)
"""
function find_logical_record(seg::SPKSegmentType2, time::Number)

    idx, Δt = divrem(time - seg.head.tstart, seg.head.tlen) 
    index = round(Int, idx)

    if index == seg.head.n 
        # This should only happen when time equals the final segment time
        index -= 1 
        Δt = seg.head.tlen 
    end 

    # Time argument for the Chebyshev polynomials
    t = Δt*seg.head.scale - 1 

    return index, t
end

"""
    get_coefficients!(daf::DAF, seg::SPKSegmentType2, desc::DAFSegmentDescriptor, index::Int)
"""
function get_coefficients!(
            daf::DAF, seg::SPKSegmentType2, desc::DAFSegmentDescriptor, index::Int)
            
    # Check whether the coefficients for this record are already loaded
    index == cache(seg).id[1] && return nothing
    cache(seg).id[1] = index 

    # Address of desired logical record 
    # (skipping mid and radius because they are all equal for SPK type 2)
    k = 8*(initial_address(desc)-1) + seg.head.recsize*index + 16

    # TODO: can we speed-up this part by casting the byte content into the array at once?
    @inbounds for j = 1:seg.head.order 
        for i = 1:seg.head.ncomp
            cache(seg).A[i, j] = get_float(
                array(daf), k + 8*(j-1) + 8*(i-1)*seg.head.order, endian(daf)
            )
        end
    end

    nothing

end

"""
    interpol(cache::SPKSegmentCache2, cheb::AbstractVector, order::Int, scale::Number)
"""
function interpol(cache::SPKSegmentCache2, cheb::AbstractVector, order::Int, scale::Number)

    len = length(cheb)

    x, y, z = 0.0, 0.0, 0.0
    @inbounds @simd for i in order+1:len
        x += cheb[i]*cache.A[1, i]
        y += cheb[i]*cache.A[2, i]
        z += cheb[i]*cache.A[3, i]
    end

    return scale*SA[x, y, z]

end

"""
    chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)
"""
function chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)

    @inbounds begin 

        cache.x1[1] = 1
        cache.x1[2] = t 
        cache.x1[3] = 2t*t - 1

        for i = 4:order 
            cache.x1[i] = 2*t*cache.x1[i-1] - cache.x1[i-2]
        end

    end

    nothing
end 

"""
    ∂chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)
"""
function ∂chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)

    @inbounds begin 

        cache.x2[1] = 0
        cache.x2[2] = 1
        cache.x2[3] = 4*t 

        for i = 4:order 
            cache.x2[i] = 2*t*cache.x2[i-1] + 2*cache.x1[i-1] - cache.x2[i-2]
        end
    end 

    nothing 
end 

"""
    ∂²chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)
"""
function ∂²chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)

    @inbounds begin 

        cache.x1[1] = 0
        cache.x1[2] = 0
        cache.x1[3] = 4 

        for i = 4:order 
            cache.x1[i] = 4*cache.x2[i-1] + 2*t*cache.x1[i-1] - cache.x2[i-2]
        end

    end

    nothing
end

"""
    ∂³chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)
"""
function ∂³chebyshev!(cache::SPKSegmentCache2, t::Number, order::Int)

    @inbounds begin 

        cache.x2[1] = 0
        cache.x2[2] = 0
        cache.x2[3] = 0 

        for i = 4:order 
            cache.x2[i] = 6*cache.x1[i-1] + 2*t*cache.x2[i-1] - cache.x2[i-2]
        end

    end

    nothing
end

