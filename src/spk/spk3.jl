
""" 
    SPKSegmentType3(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 3.
"""
function SPKSegmentType3(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader2(daf, desc)
    caches = [SPKSegmentCache2(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType3(header, caches)

end

@inline spk_field(::SPKSegmentType3) = SPK_SEGMENTLIST_MAPPING[3]


function spk_vector3(daf::DAF, seg::SPKSegmentType3, time::Number) 

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(header(seg), time)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, header(seg), cache(seg), index)

    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)
    return interpol(cache(seg), get_tmp(cache(seg).x1, time), 0, 1, 0)

end


function spk_vector6(daf::DAF, seg::SPKSegmentType3, time::Number)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(header(seg), time)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, header(seg), cache(seg), index)

    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)
    pos = interpol(cache(seg), get_tmp(cache(seg).x1, time), 0, 1, 0)
    vel = interpol(cache(seg), get_tmp(cache(seg).x1, time), 0, 1, 3)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3]
    ]

end


function spk_vector9(daf::DAF, seg::SPKSegmentType3, time::Number)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(header(seg), time)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, header(seg), cache(seg), index)

    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)
    pos = interpol(cache(seg), get_tmp(cache(seg).x1, time), 0, 1, 0)
    vel = interpol(cache(seg), get_tmp(cache(seg).x1, time), 0, 1, 3)

    # Compute 1st derivatives of Chebyshev polynomials 
    scale = seg.head.scale 
    ∂chebyshev!(cache(seg), t, seg.head.order)
    acc = interpol(cache(seg), get_tmp(cache(seg).x2, time), 1, scale, 3)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3], 
        acc[1], acc[2], acc[3]
    ]

end


function spk_vector12(daf::DAF, seg::SPKSegmentType3, time::Number)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(header(seg), time)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, header(seg), cache(seg), index)
    
    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)
    pos = interpol(cache(seg), get_tmp(cache(seg).x1, time), 0, 1, 0)
    vel = interpol(cache(seg), get_tmp(cache(seg).x1, time), 0, 1, 3)

    # Compute 1st derivatives of Chebyshev polynomials 
    scale = seg.head.scale 
    ∂chebyshev!(cache(seg), t, seg.head.order)
    acc = interpol(cache(seg), get_tmp(cache(seg).x2, time), 1, scale, 3)

    # Compute 2nd derivative of Chebyshev polynomials 
    scale *= scale
    ∂²chebyshev!(cache(seg), t, seg.head.order)
    jer = interpol(cache(seg), get_tmp(cache(seg).x1, time), 2, scale, 3)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3], 
        acc[1], acc[2], acc[3],
        jer[1], jer[2], jer[3]
    ]

end 