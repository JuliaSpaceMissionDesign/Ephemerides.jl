
""" 
    SPKSegmentHeader17(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 17.
"""
function SPKSegmentHeader17(daf::DAF, desc::DAFSegmentDescriptor)

    iaa = initial_address(desc)
    
    i0 = 8*(iaa - 1)

    # Epoch of periapsis 
    epoch = get_float(array(daf), i0, endian(daf))

    # Semi-major axis
    sma = get_float(array(daf), i0 + 8, endian(daf))

    # H term of equinoctial elements
    h = get_float(array(daf), i0 + 16, endian(daf))

    # K term of equinoctial elements
    k = get_float(array(daf), i0 + 24, endian(daf))

    # Mean longitude at epoch
    lon = get_float(array(daf), i0 + 32, endian(daf))

    # P term of equinoctial elements
    p = get_float(array(daf), i0 + 40, endian(daf))

    # Q term of equinoctial elements
    q = get_float(array(daf), i0 + 48, endian(daf))

    # Rate of longitude of periapse
    dlpdt = get_float(array(daf), i0 + 56, endian(daf))

    # Mean longitude rate (mean motion rate)
    dmldt = get_float(array(daf), i0 + 64, endian(daf))

    # Longitude of the ascending node rate
    dnodedt = get_float(array(daf), i0 + 72, endian(daf))

    # Equatorial pole right ascension
    ra = get_float(array(daf), i0 + 80, endian(daf))

    # Equatorial pole declination
    de = get_float(array(daf), i0 + 88, endian(daf))
    
    # Compute the transformation from planetary equator to the inertial reference frame
    sa, ca = sincos(ra)
    sd, cd = sincos(de)

    R = SMatrix{3, 3}(-sa, ca, 0, -ca*sd, -sa*sd, cd, ca*cd, sa*cd, sd)

    SPKSegmentHeader17(epoch, sma, h, k, lon, p, q, dlpdt, dmldt, dnodedt, ra, de, R)

end


""" 
    SPKSegmentType17(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 17.
"""
function SPKSegmentType17(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header
    header = SPKSegmentHeader17(daf, desc)
    SPKSegmentType17(header)

end

@inline spk_field(::SPKSegmentType17) = SPK_SEGMENTLIST_MAPPING[17]

function spk_vector3(::DAF, seg::SPKSegmentType17, time::Number) 

    head = header(seg)

    # Retrieve position and velocity in the orbital plane
    vf, vg, X, Y, _, _ = equinoctial_pv(head, time)

    # Compute the positions in the equinoctial reference frame 
    x = X*vf + Y*vg

    # Rotate the positions to the inertial reference frame
    return head.R*x

end

function spk_vector6(::DAF, seg::SPKSegmentType17, time::Number) 

    head = header(seg)

    # Retrieve position and velocity in the orbital plane
    vf, vg, X, Y, dX, dY = equinoctial_pv(head, time)

    # Correction factor for the periapsis rate 
    nfac = 1 - (head.dlpdt/head.dmldt)

    # Compute the rate of change of the argument of periapsis by taking the 
    # difference between the rate of the longitude of periapse and the rate of the node 
    daop = head.dlpdt - head.dnodedt

    # Correct the velocity for the precession effects
    dX = nfac*dX - daop*Y 
    dY = nfac*dY + daop*X

    # Compute the positions in the equinoctial reference frame 
    x =  X*vf +  Y*vg

    # Compute the velocity by accounting for the motion of the node 
    dnode = SA[-head.dnodedt*x[2], head.dnodedt*x[1], 0]
    v = dX*vf + dY*vg + dnode

    # Rotate the positions and velocities to the inertial reference frame
    return vcat(head.R*x, head.R*v)

end

function spk_vector9(::DAF, ::SPKSegmentType17, ::Number)
    throw(jEph.EphemerisError("Order ≥ 2 cannot be computed on segments of type 17."))
end

function spk_vector12(::DAF, ::SPKSegmentType17, ::Number)
    throw(jEph.EphemerisError("Order ≥ 2 on segments of type 17."))
end



function equinoctial_pv(head::SPKSegmentHeader17, time::Number)

    # Compute the position and velocity (the latter corrected for the perigee precession 
    # and node motion) in the orbital plane. Additionally return the coordinate axes 
    # as-well.

    # Compute the offset of the input epoch from the epoch of the elements 
    Δt = time - head.epoch

    # Compute H and K at the current epoch 
    dlp = Δt*head.dlpdt
    slp, clp = sincos(dlp)

    h =  head.h*clp + head.k*slp
    k = -head.h*slp + head.k*clp

    h2, k2 = h^2, k^2

    # Compute P and Q at the current epoch 
    node = Δt*head.dnodedt
    sn, cn = sincos(node)

    p =  head.p*cn + head.q*sn
    q = -head.p*sn + head.q*cn

    p2, q2 = p^2, q^2

    # Construct the coordinate axes 
    dI = 1/(1 + p2 + q2)

    vf = SA[1 - p2 + q2, 2*p*q, -2p]*dI 
    vg = SA[2*p*q, 1 + p2 - q2, 2q]*dI

    # Compute the mean longitude 
    λ = mod(head.lon + Δt*head.dmldt, 2π)

    # Compute the eccentric longitude from Kepler's equation
    F = kepler_eccentric_longitude(λ, h, k)
    sF, cF = sincos(F)

    # Compute Broucke's beta parameter
    b = 1/(1 + sqrt(1 - h2 - k2))

    # Compute the positions in the orbital plane 
    X = head.sma*((1 - b*h2)*cF + (h*b*sF - 1)*k)
    Y = head.sma*((1 - b*k2)*sF + (k*b*cF - 1)*h)

    # Radial distance
    r = head.sma*(1 - h*sF - k*cF)
    ra = head.dmldt*head.sma^2/r

    # Compute the velocities in the orbital plane
    dX = ra*(h*k*b*cF - (1 - h2*b)*sF)
    dY = ra*((1-k2*b)*cF - h*k*b*sF)

    return vf, vg, X, Y, dX, dY

end


function kepler_eccentric_longitude(λ, h, k)

    # Solve the equinoctial version of Kepler's equation, which is written as: 
    # λ = F + h*cos(F) - k*sin(F)

    # The algorithm is taken from: https://apps.dtic.mil/dtic/tr/fulltext/u2/a531136.pdf

    F = λ
    
    maxiter, tol = 15, 1e-16

    err = tol + 1
    iter = 0 
    while err > tol && iter < maxiter 

        sF, cF = sincos(F)
        
        fk = F + h*cF - k*sF - λ 
        dfk = 1 - h*sF - k*cF  

        ΔF = -fk/dfk 
        F += ΔF

        err = abs(ΔF)
        iter += 1
    end

    return F 

end
