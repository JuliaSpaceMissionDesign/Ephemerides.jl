
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

    # Mean longitude rate
    dmldt = get_float(array(daf), i0 + 64, endian(daf))

    # Longitude of the ascending node rate
    dnodedt = get_float(array(daf), i0 + 72, endian(daf))

    # Equatorial pole right ascension
    ra_pole = get_float(array(daf), i0 + 80, endian(daf))

    # Equatorial pole declination
    de_pole = get_float(array(daf), i0 + 88, endian(daf))
    
    SPKSegmentHeader17(epoch, sma, h, k, lon, p, q, dlpdt, dmldt, dnodedt, ra_pole, de_pole)
end

""" 
    SPKSegmentCache17(head::SPKSegmentHeader17)

Initialise the cache for an SPK segment of type 17.
"""
function SPKSegmentCache17(head::SPKSegmentHeader17) 
    SPKSegmentCache17(
        MVector(-1)
    )
end

""" 
    SPKSegmentType17(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 17.
"""
function SPKSegmentType17(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader17(daf, desc)
    caches = [SPKSegmentCache17(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType17(header, caches)

end

@inline spk_field(::SPKSegmentType17) = SPK_SEGMENTLIST_MAPPING[17]

function spk_vector3(::DAF, seg::SPKSegmentType17, time::Number) 

    head = header(seg)

    # Compute the transformation from planetary equator to the inertial reference frame
    sa, ca = sincos(head.ra_pole)
    sd, cd = sincos(head.de_pole)

    # TODO: this can be precomputed
    ROT = SMatrix{3, 3}(-sa, ca, 0, -ca*sd, -sa*sd, cd, ca*cd, sa*cd, sd)

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
    ml = head.lon + mod(Δt*head.dmldt, 2π)

    # Compute the eccentric longitude from Kepler's equation
    # F = kepler_calceph(ml, h, k)
    F = kepler_spice(ml, h, k)
    sF, cF = sincos(F)

    # Compute Broucke's beta parameter
    b = 1/(1 + sqrt(1 - h2 - k2))

    # Compute the positions in the orbital plane 
    X = head.sma*((1 - b*h2)*cF + (h*b*sF - 1)*k)
    Y = head.sma*((1 - b*k2)*sF + (k*b*cF - 1)*h)

    # Radial distance (?)
    # r = head.sma*(1 - h*sF - k*cF)

    # Compute the positions in the equinoctial reference frame 
    x = X*vf + Y*vg

    # Rotate the positions to the inertial reference frame

    # Compute the rate of change of the argument of periapse by taking the difference 
    # between the rate of the longitude of periapse and that of the node 
    # prate = head.dlpdt - head.dnodedt

    return ROT*x

end

function kepler_calceph(ml, h, k)

    # Solve the equinoctial version of Kepler's equation

    F = ml 
    
    maxiter, tol = 15, 1e-16

    err = tol + 1
    iter = 0 
    while err > tol && iter < maxiter 

        sa, ca = sincos(F)
        se = k*sa - h*ca 
        fa = F - se - ml 
        ce = k*ca + h*sa 
        f1a = 1 - ce 
        d1 = -fa/f1a 

        F = F + d1

        err = abs(d1)
        iter += 1
    end

    return F 

end

function kepler_spice(ml, h, k)

    sl, cl = sincos(ml)

    ev1 = -h*cl + k*sl 
    ev2 =  h*sl + k*cl 

    return ml + kpsolv(ev1, ev2)    

end

function kpsolv(h, k)

    y0 = -h
    xm = 0.0
    ecc = sqrt(h^2 + k^2)

    if y0 > 0 
        xu, xl = 0.0, -ecc 
    elseif y0 < 0 
        xu, xl = ecc, 0.0
    else 
        return 0.0 
    end

    maxiter = min(32, max(1, round(Int, 1/(1-ecc))))

    for _ = 1:maxiter 

        xm = max(xl, min(xu, (xl+xu)/2))
        yxm = xm - h*cos(xm) - k*sin(xm)

        if yxm > 0 
            xu = xm 
        else 
            xl = xm 
        end

    end

    # Now use Newton's method 
    x = xm 
    for _ = 1:5 

        sx, cx = sincos(x)
        yx  = x - h*cx - k*sx
        ypx = 1 + h*sx - k*cx
        x   = x - yx/ypx  

    end

    return x

end