
""" 
    SPKSegmentHeader14(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 14.
"""
function SPKSegmentHeader14(daf::DAF, desc::DAFSegmentDescriptor)

    iaa = initial_address(desc)
    faa = final_address(desc)
    
    # Read segment meta-data 
    i0 = 8*(faa-1)

    # Number of meta-data elements
    nmeta = Int(get_float(array(daf), i0, endian(daf)))

    if nmeta != 17 
        throw(ErrorException("unsupported number of SPK Type 14 meta data items."))
    end

    # Offset of the packet data from the start of a packet
    pktoff = Int(get_float(array(daf), i0 - 8, endian(daf)))

    # Size of the packets 
    pktsize = Int(get_float(array(daf), i0 - 16, endian(daf)))

    # Number of packets 
    n = Int(get_float(array(daf), i0 - 40, endian(daf)))

    # Number of directory epochs 
    ndirs = Int(get_float(array(daf), i0 - 104, endian(daf)))

    # Number of coefficients of the polynomial: this is wrongly documented in the SPK 
    # required reading documentation, which states that the polynomial degree rather than 
    # the number of coefficients is specified.
    N = get_float(array(daf), 8*(iaa-1), endian(daf))

    # Retrieve polynomial order 
    order = N - 1

    # Base address of the reference epochs  
    refbas = Int(get_float(array(daf), i0 - 88, endian(daf)))
    # Initial address of the epoch table (after all the packets)
    etid = 8*(iaa - 1 + refbas)

    # Base address of the packets (should always be 1)
    pktbas = Int(get_float(array(daf), i0 - 48, endian(daf)))
    # Initial address of the packet table (after the constants)
    ptid = 8*(iaa - 1 + pktbas)

    # Base address of the directory epochs 
    rdrbas = Int(get_float(array(daf), i0 - 112, endian(daf)))

    if ndirs == 0 
        # Load actual reference epochs 
        epochs = zeros(n)
        @inbounds for j = 1:n 
            epochs[j] = get_float(array(daf), etid + 8*(j-1), endian(daf))
        end

    else 
        # Load directory epochs 
        epochs = zeros(ndirs)
        @inbounds for j = 1:ndirs 
            epochs[j] = get_float(array(daf), 8*(iaa + rdrbas + j-2) , endian(daf))
        end
    end

    SPKSegmentHeader14(order, n, ndirs, epochs, etid, ptid, pktsize, pktoff, 6, N)
end

""" 
    SPKSegmentCache14(head::SPKSegmentHeader14)

Initialise the cache for an SPK segment of type 14.
"""
function SPKSegmentCache2(head::SPKSegmentHeader14) 
    SPKSegmentCache2(
        zeros(head.ncomp, max(3, head.N)),
        MVector(0.0, 0.0, 0.0), 
        InterpCache{Float64}(3, max(3, head.N)),
        -1
    )
end

""" 
    SPKSegmentType14(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 14.
"""
function SPKSegmentType14(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader14(daf, desc)
    caches = [SPKSegmentCache2(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType14(header, caches)

end

@inline spk_field(::SPKSegmentType14) = SPK_SEGMENTLIST_MAPPING[14]


function spk_vector3(daf::DAF, seg::SPKSegmentType14, time::Number) 

    head = header(seg)
    data = cache(seg)

    # Retrieve coefficients 
    index = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, index)

    # Normalise the time argument between [-1, 1]
    t = normalise_time(data, time)

    # Interpolate the polynomials 
    x, y, z = chebyshev(data.buff, data.A, t, 0, head.N)

    return SVector{3}(x, y, z)

end

function spk_vector6(daf::DAF, seg::SPKSegmentType14, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve coefficients 
    index = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, index)

    # Normalise the time argument between [-1, 1]
    t = normalise_time(data, time)

    # Interpolate the polynomials 
    x, y, z = chebyshev(data.buff, data.A, t, 0, head.N)
    vx, vy, vz = chebyshev(data.buff, data.A, t, 3, head.N)

    return SVector{6}(x, y, z, vx, vy, vz)

end

function spk_vector9(daf::DAF, seg::SPKSegmentType14, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve coefficients 
    index = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, index)

    # Normalise the time argument between [-1, 1]
    t = normalise_time(data, time)

    # Interpolate the polynomials 
    x, y, z = chebyshev(data.buff, data.A, t, 0, head.N)
    vx, vy, vz, ax, ay, az = ∂chebyshev(data.buff, data.A, t, 3, head.N, data.p[3])

    return SVector{9}(x, y, z, vx, vy, vz, ax, ay, az)

end

function spk_vector12(daf::DAF, seg::SPKSegmentType14, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve coefficients 
    index = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, index)
    
    # Normalise the time argument between [-1, 1]
    t = normalise_time(data, time)

    # Interpolate the polynomials 
    x, y, z = chebyshev(data.buff, data.A, t, 0, head.N)
    vx, vy, vz, ax, ay, az, jx, jy, jz = ∂²chebyshev(
        data.buff, data.A, t, 3, head.N, data.p[3]
    )

    return SVector{12}(x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz)
    
end


"""
    find_logical_record(daf::DAF, head::SPKSegmentHeader14, time::Number)
"""
function find_logical_record(daf::DAF, head::SPKSegmentHeader14, time::Number)

    # We still need to differentiate the two cases (with or without epoch directory), 
    # because in the first we are reading the final epochs from the header, in the second 
    # case, we are gradually loading them as the header contains the directory epochs! 

    # The index of the first epoch greater than time is computed in a 0-index notation 
    index = 0
    if head.ndirs == 0
        # Search through epoch table 
        while (index < head.n - 1 && time >= head.epochs[index+2])
            index += 1 
        end

    else 
        # Here, we first search through epoch directories
        subdir = 0
        @inbounds while (subdir < head.ndirs && time >= head.epochs[subdir+1])
            subdir += 1
        end

        index = max(0, 100*subdir-1)
        stop_idx = min(index + 99, head.n)

        # Find the actual epoch
        found = false
        while (index < stop_idx && !found)
            epoch = get_float(array(daf), head.etid + 8*(index+1), endian(daf))
            if time < epoch 
                found = true 
            else 
                index += 1 
            end 
        end
    end

    # The output index is returned in a 0-index notation
    return index

end

"""
    get_coefficients!(daf::DAF, head::SPKSegmentHeader14, cache::SPKSegmentCache14, index::Int)
"""
function get_coefficients!(daf::DAF, head::SPKSegmentHeader14, cache::SPKSegmentCache2, 
            index::Int)

    # Check whether the coefficients for this record are already loaded
    index == cache.id && return nothing
    cache.id = index 

    # Address of desired packet record 
    k = head.ptid + 8*head.pktoff + 8*index*(head.pktsize + head.pktoff) 

    # Retrieve the record mid point and radius 
    cache.p[1] = get_float(array(daf), k, endian(daf))
    cache.p[2] = get_float(array(daf), k + 8, endian(daf))
    cache.p[3] = 1/cache.p[2]

    k += 16 

    # # TODO: can we speed-up this part by casting the byte content into the array at once?
    @inbounds for j = 1:head.N 
        for i = 1:head.ncomp
            cache.A[i, j] = get_float(
                array(daf), k + 8*(j-1) + 8*(i-1)*head.N, endian(daf)
            )
        end
    end

    nothing

end