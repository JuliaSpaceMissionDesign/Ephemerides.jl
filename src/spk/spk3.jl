
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

function spk_vector6(daf::DAF, seg::SPKSegmentType3, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, head, data, time)

    # Normalise the time argument between [-1, 1]
    t = chebyshev_time(data, time)

    # Compute the Chebyshev polynomials
    chebyshev!(data, t, head.N)
    pos = interpol(data, get_tmp(data.x1, t), 0, 1, 0, head.N)
    vel = interpol(data, get_tmp(data.x1, t), 0, 1, 3, head.N)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3]
    ]

end


function spk_vector9(daf::DAF, seg::SPKSegmentType3, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, head, data, time)

    # Normalise the time argument between [-1, 1]
    t = chebyshev_time(data, time)

    # Compute the Chebyshev polynomials
    chebyshev!(data, t, head.N)
    pos = interpol(data, get_tmp(data.x1, t), 0, 1, 0, head.N)
    vel = interpol(data, get_tmp(data.x1, t), 0, 1, 3, head.N)

    # Compute 1st derivatives of Chebyshev polynomials 
    @inbounds scale = data.p[3]
    ∂chebyshev!(data, t, head.N)
    acc = interpol(data, get_tmp(data.x2, t), 1, scale, 3, head.N)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3], 
        acc[1], acc[2], acc[3]
    ]

end


function spk_vector12(daf::DAF, seg::SPKSegmentType3, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, head, data, time)

    # Normalise the time argument between [-1, 1]
    t = chebyshev_time(data, time)

    # Compute the Chebyshev polynomials
    chebyshev!(data, t, head.N)
    pos = interpol(data, get_tmp(data.x1, t), 0, 1, 0, head.N)
    vel = interpol(data, get_tmp(data.x1, t), 0, 1, 3, head.N)

    # Compute 1st derivatives of Chebyshev polynomials 
    @inbounds scale = data.p[3]
    ∂chebyshev!(data, t, head.N)
    acc = interpol(data, get_tmp(data.x2, t), 1, scale, 3, head.N)

    # Compute 2nd derivative of Chebyshev polynomials 
    ∂²chebyshev!(data, t, head.N)
    jer = interpol(data, get_tmp(data.x1, t), 2, scale*scale, 3, head.N)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3], 
        acc[1], acc[2], acc[3],
        jer[1], jer[2], jer[3]
    ]

end 