
""" 
    SPKSegmentHeader9(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 9 and 13.
"""
function SPKSegmentHeader9(daf::DAF, desc::DAFSegmentDescriptor)

    iaa = initial_address(desc)
    faa = final_address(desc)

    i0 = 8*(faa - 2)

    if segment_type(desc) == 9
        # Lagrange polynomials
        order = Int(get_float(array(daf), i0, endian(daf)))

        # Group size (number of states required for each interpolation)
        N = order + 1 
    else 
        # Hermite polynomials
        N = Int(get_float(array(daf), i0, endian(daf))) + 1
        order = 2N - 1
    end

    iseven = N % 2 == 0 # Check even group size 

    # Number of states stored in the record
    n = Int(get_float(array(daf), i0+8, endian(daf)))

    # Number of directory epochs (if there are 100 records, no directory epoch is stored) 
    ndirs = (n-1) ÷ 100 

    # Initial address for the epoch table (after all the states)
    etid = 8*(iaa - 1) + 48*n

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
            epochs[j] = get_float(array(daf), etid + 8*(n + j-1), endian(daf))
        end
    end
    

    SPKSegmentHeader9(n, ndirs, epochs, iaa, etid, order, N, iseven, segment_type(desc))
end

""" 
    SPKSegmentCache9(head::SPKSegmentHeader2)

Initialise the cache for an SPK segment of type 9 and 13.
"""
function SPKSegmentCache9(head::SPKSegmentHeader9) 

    if head.type == 9 
        buffsize = head.N 
        # Only for SPK 12/13 we need an array to store the the third derivatives
        nbuff = 3
    else 
        buffsize = 2*head.N
        nbuff = 4 
    end 

    SPKSegmentCache9(
        zeros(head.N),
        zeros(head.N, 6), 
        InterpCache{Float64}(nbuff, buffsize),
        MVector(-1)
    )
end

""" 
    SPKSegmentType9(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 9 and 13.
"""
function SPKSegmentType9(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader9(daf, desc)
    caches = [SPKSegmentCache9(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType9(header, caches)

end

@inline spk_field(::SPKSegmentType9) = SPK_SEGMENTLIST_MAPPING[9]


function spk_vector3(daf::DAF, seg::SPKSegmentType9, time::Number) 

    head = header(seg)
    data = cache(seg)

    # Retrieve coefficients 
    index = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, index)

    if head.type == 9
        # Interpolate the position
        x = lagrange(data.buff, data.states, data.epochs, time, 1, head.N)
        y = lagrange(data.buff, data.states, data.epochs, time, 2, head.N)
        z = lagrange(data.buff, data.states, data.epochs, time, 3, head.N)
    else 
        # Interpolate the position
        x = hermite(data.buff, data.states, data.epochs, time, 1, head.N)
        y = hermite(data.buff, data.states, data.epochs, time, 2, head.N)
        z = hermite(data.buff, data.states, data.epochs, time, 3, head.N)
    end

    return SA[x, y, z]

end

function spk_vector6(daf::DAF, seg::SPKSegmentType9, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve coefficients 
    index = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, index)

    if head.type == 9
        # Interpolate the position
        x = lagrange(data.buff, data.states, data.epochs, time, 1, head.N)
        y = lagrange(data.buff, data.states, data.epochs, time, 2, head.N)
        z = lagrange(data.buff, data.states, data.epochs, time, 3, head.N)

        # Interpolate the velocity
        vx = lagrange(data.buff, data.states, data.epochs, time, 4, head.N)
        vy = lagrange(data.buff, data.states, data.epochs, time, 5, head.N)
        vz = lagrange(data.buff, data.states, data.epochs, time, 6, head.N)
    else 
        # Interpolate the position and velocity
        x, vx = ∂hermite(data.buff, data.states, data.epochs, time, 1, head.N)
        y, vy = ∂hermite(data.buff, data.states, data.epochs, time, 2, head.N)
        z, vz = ∂hermite(data.buff, data.states, data.epochs, time, 3, head.N)
    end

    return SA[x, y, z, vx, vy, vz]

end

function spk_vector9(daf::DAF, seg::SPKSegmentType9, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve coefficients 
    index = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, index)

    if head.type == 9
        # Interpolate the position
        x = lagrange(data.buff, data.states, data.epochs, time, 1, head.N)
        y = lagrange(data.buff, data.states, data.epochs, time, 2, head.N)
        z = lagrange(data.buff, data.states, data.epochs, time, 3, head.N)

        # Interpolate the velocity and acceleration
        vx, ax = ∂lagrange(data.buff, data.states, data.epochs, time, 4, head.N)
        vy, ay = ∂lagrange(data.buff, data.states, data.epochs, time, 5, head.N)
        vz, az = ∂lagrange(data.buff, data.states, data.epochs, time, 6, head.N)
    else 
        # Interpolate the position, velocity and acceleration
        x, vx, ax = ∂²hermite(data.buff, data.states, data.epochs, time, 1, head.N)
        y, vy, ay = ∂²hermite(data.buff, data.states, data.epochs, time, 2, head.N)
        z, vz, az = ∂²hermite(data.buff, data.states, data.epochs, time, 3, head.N)
    end

    return SA[x, y, z, vx, vy, vz, ax, ay, az]

end

function spk_vector12(daf::DAF, seg::SPKSegmentType9, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve coefficients 
    index = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, index)

    if head.type == 9
        # Interpolate the position
        x = lagrange(data.buff, data.states, data.epochs, time, 1, head.N)
        y = lagrange(data.buff, data.states, data.epochs, time, 2, head.N)
        z = lagrange(data.buff, data.states, data.epochs, time, 3, head.N)
        
        # Interpolate the velocity, acceleration and jerk
        vx, ax, jx = ∂²lagrange(data.buff, data.states, data.epochs, time, 4, head.N)
        vy, ay, jy = ∂²lagrange(data.buff, data.states, data.epochs, time, 5, head.N)
        vz, az, jz = ∂²lagrange(data.buff, data.states, data.epochs, time, 6, head.N)
    else 
        # Interpolate the position, velocity, acceleration and jerk
        x, vx, ax, jx = ∂³hermite(data.buff, data.states, data.epochs, time, 1, head.N)
        y, vy, ay, jy = ∂³hermite(data.buff, data.states, data.epochs, time, 2, head.N)
        z, vz, az, jz = ∂³hermite(data.buff, data.states, data.epochs, time, 3, head.N)
    end

    return SA[x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz]
    
end


"""
    find_logical_record(daf::DAF, head::SPKSegmentHeader1, time::Number)
"""
function find_logical_record(daf::DAF, head::SPKSegmentHeader9, time::Number)

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
        stop_idx = min(index + 100, head.n)

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

    if head.iseven 
        first = index - (head.N ÷ 2)
    elseif index == 0 
        first = 0
    else
        # Need to check which one of the two epochs is closer 
        if head.ndirs == 0 
            @inbounds e1, e2 = head.epochs[index], head.epochs[index + 1]
        else 
            e1 = get_float(array(daf), head.etid + 8*(index-1), endian(daf))
            e2 = get_float(array(daf), head.etid + 8*index, endian(daf))
        end
            
        if time - e1 < e2 - time
            index = index - 1              
        end

        first = index - (head.N - 1) ÷ 2
    end

    # The output index is returned in a 0-index notation
    return min(max(0, first), head.n - head.N)

end

"""
    get_coefficients!(daf::DAF, head::SPKSegmentHeader9, cache::SPKSegmentCache9, index::Int)
"""
function get_coefficients!(daf::DAF, head::SPKSegmentHeader9, cache::SPKSegmentCache9, 
            index::Int)

    # Check whether the coefficients for this record are already loaded
    index == cache.id[1] && return nothing
    cache.id[1] = index 

    # Address of desired logical record 
    i0 = 8*(head.iaa - 1) + 48*index # 48 = 6 (states) * 8 bytes

    # TODO: can we speed-up this part by casting the byte content into the array at once?
    
    # Here we load-up the state coefficients
    @inbounds for j = 1:6
        for i = 1:head.N 
            cache.states[i, j] = get_float(array(daf), i0 + 48*(i-1) + 8*(j-1), endian(daf))
        end
    end

    # Here we load the epoch values 
    if head.ndirs == 0 
        # Extract the epochs directly from the ones stored in the header since there is 
        # no directory epoch
        @inbounds for i = 1:head.N 
            cache.epochs[i] = head.epochs[index+i]
        end

    else 
        # Read the epochs from the binary file 
        i0 = 8*(head.iaa - 1) + 48*head.n + 8*index
        @inbounds for i = 1:head.N 
            cache.epochs[i] = get_float(array(daf), i0 + 8*(i-1), endian(daf))
        end
    end

    nothing 
end