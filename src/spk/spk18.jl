
""" 
    SPKSegmentHeader18(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 18.
"""
function SPKSegmentHeader18(daf::DAF, desc::DAFSegmentDescriptor)

    iaa = initial_address(desc)
    faa = final_address(desc)

    i0 = 8*(faa - 3)

    # Retrieve the segment subtype, currently only 0 (Hermite) or 1 (Lagrange) exist
    subtype = Int(get_float(array(daf), i0, endian(daf)))

    # Retrieve window size (i.e., number of interpolating points)
    N = Int(get_float(array(daf), i0 + 8, endian(daf)))

    if isodd(N) 
        throw(jEph.EphemerisError(
            "Unexpected odd window size $N found in SPK type 18 (must be even).")
        )
    end

    # Retrieve number of packets 
    n = Int(get_float(array(daf), i0 + 16, endian(daf)))

    if n < 2
        throw(jEph.EphemerisError(
            "The packet count in the SPK Type 18 is below the minimum allowed value (2)")
        )
    end

    # Compute polynomial order and packet size
    if subtype == 0 # Hermite polynomial 
        order = 2N - 1
        packetsize = 12 
    elseif subtype == 1 # Lagrange polynomial
        order = N - 1
        packetsize = 6
    else 
        throw(jEph.EphemerisError("Unexpected subtype $subtype found in SPK type 18"))
    end 

    # Number of directory epochs (if there are 100 packets, no directory epoch is stored) 
    ndirs = (n-1) ÷ 100

    # Initial address for the epoch table (after all the states)
    etid = 8*(iaa - 1) + 8*packetsize*n
    
    if ndirs == 0 
        # Load actual epochs 
        epochs = zeros(n)
        @inbounds for j = 1:n 
            epochs[j] = get_float(array(daf), etid + 8*(j-1), endian(daf))
        end
    else 
        # Load directory epochs 
        epochs = zeros(ndirs)
        @inbounds for j = 1:ndirs
            epochs[j] = get_float(array(daf), etid + 8*(n+j-1), endian(daf))
        end
    end

    SPKSegmentHeader18(n, ndirs, epochs, iaa, etid, order, N, subtype, packetsize)
end

""" 
    SPKSegmentCache18(spkhead::SPKSegmentHeader2)

Initialise the cache for an SPK segment of type 18.
"""
function SPKSegmentCache18(header::SPKSegmentHeader18) 
    
    if header.subtype == 0 # Hermite
        buffsize = 2*header.N  
    else
        buffsize = header.N # Lagrange
    end

    # For this spk type, both subtypes share the same number of buffers because 
    # also for the Hermite polynomials we have two independent set of coefficients 
    # for the position and the velocity
    nbuff = 3
    
    SPKSegmentCache18(
        MVector(-1, -1, header.N),
        zeros(header.N), 
        zeros(header.N, header.packetsize),
        InterpCache{Float64}(nbuff, buffsize)
    )
end

""" 
    SPKSegmentType18(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 18.
"""
function SPKSegmentType18(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader18(daf, desc)
    caches = [SPKSegmentCache18(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType18(header, caches)

end

@inline spk_field(::SPKSegmentType18) = SPK_SEGMENTLIST_MAPPING[18]


function spk_vector3(daf::DAF, seg::SPKSegmentType18, time::Number) 

    head = header(seg)
    data = cache(seg)

    # Retrieve the indexes of the first and last point and the size of the 
    # polynomial sliding-window. 
    first, last = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, first, last)

    # Retrieve the windowsize
    @inbounds N = data.p[3]

    # Interpolate the position coefficients
    if head.subtype == 0
        x = hermite(data.buff, data.states, data.epochs, time, 1, N)
        y = hermite(data.buff, data.states, data.epochs, time, 2, N)
        z = hermite(data.buff, data.states, data.epochs, time, 3, N)
    else 
        x = lagrange(data.buff, data.states, data.epochs, time, 1, N)
        y = lagrange(data.buff, data.states, data.epochs, time, 2, N)
        z = lagrange(data.buff, data.states, data.epochs, time, 3, N)
    end

    return SA[x, y, z]

end

function spk_vector6(daf::DAF, seg::SPKSegmentType18, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve the indexes of the first and last point and the size of the 
    # polynomial sliding-window. 
    first, last = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, first, last)

    # Retrieve the windowsize
    @inbounds N = data.p[3]

    if head.subtype == 0
        # Interpolate the position
        x = hermite(data.buff, data.states, data.epochs, time, 1, N)
        y = hermite(data.buff, data.states, data.epochs, time, 2, N)
        z = hermite(data.buff, data.states, data.epochs, time, 3, N)

        # Interpolate the velocity
        vx = hermite(data.buff, data.states, data.epochs, time, 7, N)
        vy = hermite(data.buff, data.states, data.epochs, time, 8, N)
        vz = hermite(data.buff, data.states, data.epochs, time, 9, N)

    else 
        # Interpolate the position
        x = lagrange(data.buff, data.states, data.epochs, time, 1, N)
        y = lagrange(data.buff, data.states, data.epochs, time, 2, N)
        z = lagrange(data.buff, data.states, data.epochs, time, 3, N)

        # Interpolate the velocity
        vx = lagrange(data.buff, data.states, data.epochs, time, 4, N)
        vy = lagrange(data.buff, data.states, data.epochs, time, 5, N)
        vz = lagrange(data.buff, data.states, data.epochs, time, 6, N)
    end

    return SA[x, y, z, vx, vy, vz]

end

function spk_vector9(daf::DAF, seg::SPKSegmentType18, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve the indexes of the first and last point and the size of the 
    # polynomial sliding-window. 
    first, last = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, first, last)

    # Retrieve the windowsize
    @inbounds N = data.p[3]

    if head.subtype == 0
        # Interpolate the position
        x = hermite(data.buff, data.states, data.epochs, time, 1, N)
        y = hermite(data.buff, data.states, data.epochs, time, 2, N)
        z = hermite(data.buff, data.states, data.epochs, time, 3, N)

        # Interpolate the velocity and acceleration
        vx, ax = ∂hermite(data.buff, data.states, data.epochs, time, 7, N)
        vy, ay = ∂hermite(data.buff, data.states, data.epochs, time, 8, N)
        vz, az = ∂hermite(data.buff, data.states, data.epochs, time, 9, N)

    else 
        # Interpolate the position
        x = lagrange(data.buff, data.states, data.epochs, time, 1, N)
        y = lagrange(data.buff, data.states, data.epochs, time, 2, N)
        z = lagrange(data.buff, data.states, data.epochs, time, 3, N)

        # Interpolate the velocity and acceleration
        vx, ax = ∂lagrange(data.buff, data.states, data.epochs, time, 4, N)
        vy, ay = ∂lagrange(data.buff, data.states, data.epochs, time, 5, N)
        vz, az = ∂lagrange(data.buff, data.states, data.epochs, time, 6, N)
    end

    return SA[x, y, z, vx, vy, vz, ax, ay, az]

end

function spk_vector12(daf::DAF, seg::SPKSegmentType18, time::Number)
    head = header(seg)
    data = cache(seg)

    # Retrieve the indexes of the first and last point and the size of the 
    # polynomial sliding-window. 
    first, last = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, first, last)

    # Retrieve the windowsize
    @inbounds N = data.p[3]

    if head.subtype == 0
        # Interpolate the position
        x = hermite(data.buff, data.states, data.epochs, time, 1, N)
        y = hermite(data.buff, data.states, data.epochs, time, 2, N)
        z = hermite(data.buff, data.states, data.epochs, time, 3, N)

        # Interpolate the velocity, acceleration and jerk
        vx, ax, jx = ∂²hermite(data.buff, data.states, data.epochs, time, 7, N)
        vy, ay, jy = ∂²hermite(data.buff, data.states, data.epochs, time, 8, N)
        vz, az, jz = ∂²hermite(data.buff, data.states, data.epochs, time, 9, N)

    else 
        # Interpolate the position
        x = lagrange(data.buff, data.states, data.epochs, time, 1, N)
        y = lagrange(data.buff, data.states, data.epochs, time, 2, N)
        z = lagrange(data.buff, data.states, data.epochs, time, 3, N)

        # Interpolate the velocity, acceleration and jerk
        vx, ax, jx = ∂²lagrange(data.buff, data.states, data.epochs, time, 4, N)
        vy, ay, jy = ∂²lagrange(data.buff, data.states, data.epochs, time, 5, N)
        vz, az, jz = ∂²lagrange(data.buff, data.states, data.epochs, time, 6, N)
    end

    return SA[x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz]

end


"""
    find_logical_record(daf::DAF, head::SPKSegmentHeader18, time::Number)
"""
function find_logical_record(daf::DAF, head::SPKSegmentHeader18, time::Number)

    # We still need to differentiate the two cases (with or without epoch directory), 
    # because in the first we are reading the final epochs from the header, in the second 
    # case, we are gradually loading them as the header contains the directory epochs! 

    # The index of the first epoch greater than time is computed in a 0-index notation 
    index = 0
    if head.ndirs == 0
        # Search through epoch table 
        @inbounds while (index < head.n - 1 && head.epochs[index+1] < time)
            index += 1 
        end
    else 
        # Here, we first search through epoch directories
        subdir = 0
        @inbounds while (subdir < head.ndirs && head.epochs[subdir+1] < time)
            subdir += 1
        end

        index = subdir*100
        stop_idx = max(index + 100, head.n)

        # Find the actual epoch
        found = false
        while (index < stop_idx - 1 && !found)
            epoch = get_float(array(daf), head.etid + 8*index, endian(daf))
            
            if epoch >= time 
                found = true 
            else 
                index += 1 
            end 
        end
    end

    # At this point of the code we have found the index of the first epoch greater than 
    # the requested time.

    # Half window size (by design N is even)
    hsize = head.N ÷ 2

    # Compute lower and upper indexes (in 0-based notation) of the segment points,
    # differently from types 8/9/12/13 here the window size is variable when the time 
    # is close to the segment boundaries 
    first = max(0, index - hsize)
    last = min(head.n - 1, index + hsize - 1)

    return first, last

end


function get_coefficients!(daf::DAF, head::SPKSegmentHeader18, cache::SPKSegmentCache18, 
            first::Int, last::Int)

    # Since the window might be truncated at one side, we need to check for both 
    # boundaries to see whether the coefficients have already been loaded 

    # Check whether the packets for this polynomial are already loaded 
    if first == cache.p[1] && last == cache.p[2]
        return nothing 
    end 

    @inbounds cache.p[1] = first 
    @inbounds cache.p[2] = last 

    # Compute the actual window size
    @inbounds cache.p[3] = last - first + 1

    # Address of the first desired logical record 
    i0 = 8*(head.iaa - 1) + 8*head.packetsize*first

    # Here we load-up the state coefficients 
    @inbounds for j = 1:head.packetsize
        for i = 1:cache.p[3]
            cache.states[i, j] = get_float(
                array(daf), i0 + 8*head.packetsize*(i-1) + 8*(j-1), endian(daf)
            )
        end
    end

    # Here we load the epoch values 
    if head.ndirs == 0 
        # Extract the epochs directly from the ones stored in the header since there 
        # is no directory epoch 
        @inbounds for i = 1:cache.p[3]
            cache.epochs[i] = head.epochs[first+i]
        end
    else 
        # Read the epochs from the binary file
        i0 = 8*(head.iaa - 1) + 8*head.packetsize*head.n + 8*first 
        @inbounds for i = 1:cache.p[3]
            cache.epochs[i] = get_float(array(daf), i0 + 8*(i-1), endian(daf))
        end
    end

    nothing

end 
