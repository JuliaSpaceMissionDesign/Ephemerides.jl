
""" 
    SPKSegmentHeader8(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 8 and 12.
"""
function SPKSegmentHeader8(daf::DAF, desc::DAFSegmentDescriptor)

    iaa = initial_address(desc)
    faa = final_address(desc)

    i0 = 8*(faa - 4)

    # Retrieve segment starting epoch and length 
    tstart = get_float(array(daf), i0, endian(daf))
    tlen = get_float(array(daf), i0+8, endian(daf))

    # Retrieve polynomial order 
    if segment_type(desc) == 8
        # Lagrange polynomials
        order = Int(get_float(array(daf), i0 + 16, endian(daf)))
        N = order + 1
    else
        # Hermite polynomials
        N = Int(get_float(array(daf), i0 + 16, endian(daf))) + 1
        order = 2N - 1
    end 

    # Number of states stored in the record 
    n = Int(get_float(array(daf), i0 + 24, endian(daf)))

    # Check even group size 
    iseven = N % 2 == 0

    SPKSegmentHeader8(tstart, tlen, order, N, n, iaa, iseven, segment_type(desc))
end

""" 
    SPKSegmentCache8(head::SPKSegmentHeader2)

Initialise the cache for an SPK segment of type 8 and 12.
"""
function SPKSegmentCache8(head::SPKSegmentHeader8) 

    if head.type == 8 
        buffsize = head.N 
        nbuff = 3
    else 
        buffsize = 2*head.N
        nbuff = 4
    end 

    SPKSegmentCache8(
        zeros(head.N, 6), 
        InterpCache{Float64}(nbuff, buffsize),
        MVector(-1)
    )
end

""" 
    SPKSegmentType8(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 8 and 12.
"""
function SPKSegmentType8(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader8(daf, desc)
    caches = [SPKSegmentCache8(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType8(header, caches)

end

@inline spk_field(::SPKSegmentType8) = SPK_SEGMENTLIST_MAPPING[8]


function spk_vector3(daf::DAF, seg::SPKSegmentType8, time::Number) 

    head = header(seg)
    data = cache(seg)

    # Retrieve Lagrange coefficients 
    index = find_logical_record(head, time)
    get_coefficients!(daf, head, data, index)

    # Normalise the time argument between [1, 2]
    Δt = normalise_time(head, time, index)

    # Interpolate the position
    if head.type == 8
        x = lagrange(data.buff, data.states, Δt, 1, head.N)
        y = lagrange(data.buff, data.states, Δt, 2, head.N)
        z = lagrange(data.buff, data.states, Δt, 3, head.N)
    else 
        x = hermite(data.buff, data.states, Δt, 1, head.N, head.tlen)
        y = hermite(data.buff, data.states, Δt, 2, head.N, head.tlen)
        z = hermite(data.buff, data.states, Δt, 3, head.N, head.tlen)
    end

    return SA[x, y, z]

end

function spk_vector6(daf::DAF, seg::SPKSegmentType8, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve Lagrange coefficients 
    index = find_logical_record(head, time)
    get_coefficients!(daf, head, data, index)

    # Normalise the time argument between [1, 2]
    Δt = normalise_time(head, time, index)

    if head.type == 8
        # Interpolate the position
        x = lagrange(data.buff, data.states, Δt, 1, head.N)
        y = lagrange(data.buff, data.states, Δt, 2, head.N)
        z = lagrange(data.buff, data.states, Δt, 3, head.N)

        # Interpolate the velocity
        vx = lagrange(data.buff, data.states, Δt, 4, head.N)
        vy = lagrange(data.buff, data.states, Δt, 5, head.N)
        vz = lagrange(data.buff, data.states, Δt, 6, head.N)
    else 
        # Interpolate the position and velocity
        x, vx = ∂hermite(data.buff, data.states, Δt, 1, head.N, head.tlen)
        y, vy = ∂hermite(data.buff, data.states, Δt, 2, head.N, head.tlen)
        z, vz = ∂hermite(data.buff, data.states, Δt, 3, head.N, head.tlen)
    end

    return SA[x, y, z, vx, vy, vz]

end

function spk_vector9(daf::DAF, seg::SPKSegmentType8, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve Lagrange coefficients 
    index = find_logical_record(head, time)
    get_coefficients!(daf, head, data, index)

    # Normalise the time argument between [1, 2]
    Δt = normalise_time(head, time, index)

    if head.type == 8
        # Interpolate the position
        x = lagrange(data.buff, data.states, Δt, 1, head.N)
        y = lagrange(data.buff, data.states, Δt, 2, head.N)
        z = lagrange(data.buff, data.states, Δt, 3, head.N)
        
        # Interpolate the velocity and acceleration
        vx, ax = ∂lagrange(data.buff, data.states, Δt, 4, head.N, head.tlen)
        vy, ay = ∂lagrange(data.buff, data.states, Δt, 5, head.N, head.tlen)
        vz, az = ∂lagrange(data.buff, data.states, Δt, 6, head.N, head.tlen)
    else 
        # Interpolate the position, velocity and acceleration
        x, vx, ax = ∂²hermite(data.buff, data.states, Δt, 1, head.N, head.tlen)
        y, vy, ay = ∂²hermite(data.buff, data.states, Δt, 2, head.N, head.tlen)
        z, vz, az = ∂²hermite(data.buff, data.states, Δt, 3, head.N, head.tlen)
    end

    return SA[x, y, z, vx, vy, vz, ax, ay, az]

end

function spk_vector12(daf::DAF, seg::SPKSegmentType8, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve Lagrange coefficients 
    index = find_logical_record(head, time)
    get_coefficients!(daf, head, data, index)

    # Normalise the time argument between [1, 2]
    Δt = normalise_time(head, time, index)

    if head.type == 8
        # Interpolate the position
        x = lagrange(data.buff, data.states, Δt, 1, head.N)
        y = lagrange(data.buff, data.states, Δt, 2, head.N)
        z = lagrange(data.buff, data.states, Δt, 3, head.N)
        
        # Interpolate the velocity, acceleration and jerk
        vx, ax, jx = ∂²lagrange(data.buff, data.states, Δt, 4, head.N, head.tlen)
        vy, ay, jy = ∂²lagrange(data.buff, data.states, Δt, 5, head.N, head.tlen)
        vz, az, jz = ∂²lagrange(data.buff, data.states, Δt, 6, head.N, head.tlen)
    else 
        # Interpolate the position, velocity, acceleration and jerk
        x, vx, ax, jx = ∂³hermite(data.buff, data.states, Δt, 1, head.N, head.tlen)
        y, vy, ay, jy = ∂³hermite(data.buff, data.states, Δt, 2, head.N, head.tlen)
        z, vz, az, jz = ∂³hermite(data.buff, data.states, Δt, 3, head.N, head.tlen)
    end
    
    return SA[x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz]
end


"""
    find_logical_record(head::SPKSegmentHeader8, time::Number)
"""
function find_logical_record(head::SPKSegmentHeader8, time::Number)

    Δt = time - head.tstart 

    if head.iseven # even group size 
        low = Int(Δt ÷ head.tlen) + 1
        first = low - (head.N ÷ 2) + 1
    else # odd group size  
        near = round(Int, Δt/head.tlen) + 1
        first = near - (head.N - 1) ÷ 2
    end

    index = min(max(1, first), head.n - head.N + 1)
    return index 

end

"""
    get_coefficients!(daf::DAF, head::SPKSegmentHeader8, cache::SPKSegmentCache8, index::Int)
"""
function get_coefficients!(daf::DAF, head::SPKSegmentHeader8, cache::SPKSegmentCache8, 
            index::Int)

    # Check whether the coefficients for this record are already loaded
    index == cache.id[1] && return nothing
    cache.id[1] = index 

    # Address of desired logical record 
    i0 = 8*(head.iaa - 1) + 48*(index-1) # 6*8*(index - 1)

    # TODO: can we speed-up this part by casting the byte content into the array at once?
    @inbounds for j = 1:6
        for i = 1:head.N 
            cache.states[i, j] = get_float(
                array(daf), i0 + 48*(i-1) + 8*(j-1), endian(daf)
            )
        end
    end

end

"""
    normalise_time(head::SPKSegmentHeader8, time::Number, index::Int)

Returned a normalised time that starts at 1 at the beginning of the interval.
"""
function normalise_time(head::SPKSegmentHeader8, time::Number, index::Int)
    tbegin = head.tstart + (index-1)*head.tlen 
    return (time - tbegin)/head.tlen + 1
end

