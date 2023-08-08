
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
    pos = _inner_vector3(seg, t)
    vel = interpol(cache(seg), cache(seg).x1, 0, 1, 3)

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
    vel = interpol(cache(seg), cache(seg).x1, 0, 1, 3)

    # Compute 1st derivatives of Chebyshev polynomials 
    scale = seg.head.scale 
    ∂chebyshev!(cache(seg), t, seg.head.order)
    acc = interpol(cache(seg), cache(seg).x2, 1, scale, 3)

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
    vel = interpol(cache(seg), cache(seg).x1, 0, 1, 3)

    # Compute 1st derivatives of Chebyshev polynomials 
    scale = seg.head.scale 
    ∂chebyshev!(cache(seg), t, seg.head.order)
    acc = interpol(cache(seg), cache(seg).x2, 1, scale, 3)

    # Compute 2nd derivative of Chebyshev polynomials 
    scale *= scale
    ∂²chebyshev!(cache(seg), t, seg.head.order)
    jer = interpol(cache(seg), cache(seg).x1, 2, scale, 3)

    return @inbounds SA[
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3], 
        acc[1], acc[2], acc[3],
        jer[1], jer[2], jer[3]
    ]

end

function _inner_vector3(seg::SPKSegmentType3, ::Number)

    # Interpolate the body position
    return interpol(cache(seg), cache(seg).x1, 0, 1, 0)

end

# --------------------
#  ForwardDiff Rules
# --------------------

function _inner_vector3(seg::SPKSegmentType3, t::Dual{T}) where T 

    pos = _inner_vector3(seg, value(t))
    vel = interpol(cache(seg), cache(seg).x1, 0, 1, 3)

    return Dual{T}(pos, vel*partials(time))

end

