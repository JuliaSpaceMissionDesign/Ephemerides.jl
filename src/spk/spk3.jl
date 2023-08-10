
""" 
    SPKSegmentType3(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 3.
"""
function SPKSegmentType3(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    head = SPKSegmentHeader2(daf, desc)
    cache = SPKSegmentCache2(head)

    SPKSegmentType3(head, cache)

end

@inline spk_field(::SPKSegmentType3) = SPK_SEGMENTLIST_MAPPING[3]


function spk_vector3(daf::DAF, seg::SPKSegmentType3, desc::DAFSegmentDescriptor, time::Number) 

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(header(seg), time)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, header(seg), cache(seg), desc, index)

    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)
    return _inner_vector3(seg, t)

end


function spk_vector6(daf::DAF, seg::SPKSegmentType3, desc::DAFSegmentDescriptor, time::Number)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(header(seg), time)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, header(seg), cache(seg), desc, index)

    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)
    # pos = _inner_vector3(seg, t)
    pos = interpol(cache(seg), get_tmp(cache(seg).x1, time), 0, 1, 0)
    vel = interpol(cache(seg), get_tmp(cache(seg).x1, time), 0, 1, 3)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3]
    ]

end

function spk_vector9(daf::DAF, seg::SPKSegmentType3, desc::DAFSegmentDescriptor, time::Number)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(header(seg), time)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, header(seg), cache(seg), desc, index)

    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)
    pos = _inner_vector3(seg, t)
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

function spk_vector12(daf::DAF, seg::SPKSegmentType3, desc::DAFSegmentDescriptor, time::Number)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(header(seg), time)

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, header(seg), cache(seg), desc, index)
    
    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)
    pos = _inner_vector3(seg, t)
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

function _inner_vector3(seg::SPKSegmentType3, t::Number)

    # Interpolate the body position
    return interpol(cache(seg), get_tmp(cache(seg).x1, t), 0, 1, 0)

end

# --------------------
#  ForwardDiff Rules
# --------------------

function _inner_vector3(seg::SPKSegmentType3, t::Dual{T}) where T 

    pos = interpol(cache(seg), get_tmp(cache(seg).x1, value(t)), 0, 1, 0)
    vel = interpol(cache(seg), get_tmp(cache(seg).x1, value(t)), 0, 1, 3)
    
    println(pos, "\n", vel)

    # return vel

    return SA[
        Dual{T}(pos[1], vel[1]*partials(t)),
        Dual{T}(pos[2], vel[2]*partials(t)), 
        Dual{T}(pos[3], vel[3]*partials(t)), 
    ] 

end


function spk_vector4(daf::DAF, seg::SPKSegmentType3, desc::DAFSegmentDescriptor, time::Dual{T}) where {T} 

    # Find the logical record containing the Chebyshev coefficients at `time`
    index, t = find_logical_record(header(seg), value(time))

    # Retrieve Chebyshev coefficients 
    get_coefficients!(daf, header(seg), cache(seg), desc, index)

    # Compute the Chebyshev polynomials
    chebyshev!(cache(seg), t, seg.head.order)

    pos = _inner_vector3(seg, t)
    vel = interpol(cache(seg), get_tmp(cache(seg).x1, t), 0, 1, 3)
    
    return SA[
        Dual{T}(pos[1], vel[1]*partials(time)), 
        Dual{T}(pos[2], vel[2]*partials(time)), 
        Dual{T}(pos[3], vel[3]*partials(time)), 
    ]

end