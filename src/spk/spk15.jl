
""" 
    SPKSegmentHeader15(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 15.
"""
function SPKSegmentHeader15(daf::DAF, desc::DAFSegmentDescriptor)

    iaa = initial_address(desc)
    
    # Read segment data 
    i0 = 8*(iaa-1)

    # Epoch of periapsis
    epoch = get_float(array(daf), i0, endian(daf))

    # Trajectory pole
    tp = zeros(3) 
    @inbounds for j = 1:3 
        tp[j] = get_float(array(daf), i0 + 8j, endian(daf))
    end
    
    # Normalise the vector 
    tp[:] = vhat(tp)
    i0 += 24

    # Periapsis vector
    pa = zeros(3)
    @inbounds for j = 1:3 
        pa[j] = get_float(array(daf), i0 + 8j, endian(daf))
    end
    
    # Normalise the vector 
    pa[:] = vhat(pa)
    i0 += 24 

    # Semi-latus rectum
    p = get_float(array(daf), i0 + 8, endian(daf))

    # Eccentricity 
    ecc = get_float(array(daf), i0 + 16, endian(daf))

    # J2 processing flag 
    j2f = get_float(array(daf), i0 + 24, endian(daf))
    i0 += 24
    
    # Central Body pole 
    pv = zeros(3)
    @inbounds for j = 1:3 
        pv[j] = get_float(array(daf), i0 + 8j, endian(daf))
    end

    # Normalise the vector 
    pv[:] = vhat(pv)
    i0 += 24 

    # Central Body Gravitational parameter 
    GM = get_float(array(daf), i0 + 8, endian(daf))

    # Central Body J2 
    J2 = get_float(array(daf), i0 + 16, endian(daf))

    # Central Body equatorial radius 
    R = get_float(array(daf), i0 + 24, endian(daf))

    # Even if the J2 perturbations are requested, check whether the orbit satisfies the 
    # requirements for it to be applied (i.e., is elliptic and does not impact the planet)
    vj2 = j2f != 3 && ecc < 1 && p/(1 + ecc) > R && J2 != 0 

    if vj2 
        # Compute the mean anomaly rate of change 
        x = 1 - ecc^2 
        dmdt = (x/p)*sqrt(GM*x/p)

        # Compute cos(inc)
        cosinc = vdot(pv, tp)

        # Compute gain factors for regression of node and precession of pericenter
        kz = 3/2*J2*(R/p)^2
        kn = -kz*cosinc
        kp = kz*(5/2*cosinc^2 - 1/2)
    else 
        dmdt = 0.0 
        kn, kp = 0.0, 0.0
    end

    SPKSegmentHeader15(
        epoch, SVector{3}(tp), SVector{3}(pv), SVector{3}(pa),
        p, ecc, j2f, vj2, GM, J2, R, dmdt, kn, kp
    )

end

""" 
    SPKSegmentType15(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 15.
"""
function SPKSegmentType15(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header
    head = SPKSegmentHeader15(daf, desc)

    # Compute the distance and the speed at periapse 
    rp = head.p/(1 + head.ecc)
    vp = sqrt(head.GM/head.p)*(1 + head.ecc)

    # Get the position and velocity vectors at the periapse 
    pos = rp*head.pa 
    vel = vp*vcross(head.tp, head.pa)

    # Initialise the two-body caches
    caches = [TwoBodyUniversalCache(head.GM, pos, vel) for _ in 1:Threads.nthreads()]
    SPKSegmentType15(head, caches)

end

@inline spk_field(::SPKSegmentType15) = SPK_SEGMENTLIST_MAPPING[15]

function spk_vector3(::DAF, seg::SPKSegmentType15, time::Number) 

    head = header(seg)
    data = cache(seg)

    # Propagate the state at periapses to the requested epoch 
    Δt = time - head.epoch 

    pos, _ = propagate_twobody(data, Δt)
    if head.vj2  

        # Account for the J₂ effect. The rotations can be performed in any order, but if 
        # the line of node regressions was accounted frist, we would also have to rotate 
        # the pole before precessing the line of apsides.

        # Determine how far the line of nodes and periapsis have moved 
        ΔΩ, Δω = compute_j2_effects(head, pos, Δt)

        # Precesses the periapsis by rotating the state vector about the trajectory pole 
        if head.j2f != 1 
            pos = vrot(pos, head.tp, Δω)
        end

        # Regresses the line of nodes by rotating the state about the pole of the body 
        if head.j2f != 2 
            pos = vrot(pos, head.pv, ΔΩ)
        end

    end

    return pos

end

function spk_vector6(::DAF, seg::SPKSegmentType15, time::Number) 

    head = header(seg)
    data = cache(seg)

    # Propagate the state at periapses to the requested epoch 
    Δt = time - head.epoch 

    pos, vel = propagate_twobody(data, Δt)
    if head.vj2  

        # Account for the J₂ effect. The rotations can be performed in any order, but if 
        # the line of node regressions was accounted frist, we would also have to rotate 
        # the pole before precessing the line of apsides.

        # Determine how far the line of nodes and periapsis have moved 
        ΔΩ, Δω = compute_j2_effects(head, pos, Δt)

        # Precesses the periapsis by rotating the state vector about the trajectory pole 
        if head.j2f != 1 
            pos = vrot(pos, head.tp, Δω)
            vel = vrot(vel, head.tp, Δω)
        end

        # Regresses the line of nodes by rotating the state about the pole of the body 
        if head.j2f != 2 
            pos = vrot(pos, head.pv, ΔΩ)
            vel = vrot(vel, head.pv, ΔΩ)
        end

    end

    return vcat(pos, vel)

end


function spk_vector9(::DAF, ::SPKSegmentType15, ::Number)
    throw(jEph.EphemerisError("Order ≥ 2 cannot be computed on segments of type 15."))
end

function spk_vector12(::DAF, ::SPKSegmentType15, ::Number)
    throw(jEph.EphemerisError("Order ≥ 2 on segments of type 15."))
end


function compute_j2_effects(head::SPKSegmentHeader15, pos, Δt::Number)

    # Compute the change in mean anomaly since periapsis 
    ΔM = head.dmdt*Δt

    # Compute the angle θ such that -π < θ < π and ΔM = θ + 2kπ 
    θ = mod(ΔM, 2π)
    if abs(θ) > π 
        θ -= 2π*sign(θ)
    end 

    # Get the accumulated true anomaly from the propagated state θ and the accumulated 
    # mean anomaly prior to this orbit 
    TA = vsep(head.pa, pos)*sign(θ) + ΔM - θ

    # Determine how far the line of nodes and periapsis have moved 
    ΔΩ = TA*head.kn 
    Δω = TA*head.kp

    return ΔΩ, Δω

end