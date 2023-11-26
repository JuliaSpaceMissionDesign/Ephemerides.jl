
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

