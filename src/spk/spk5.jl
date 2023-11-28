
""" 
    SPKSegmentHeader5(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 5.
"""
function SPKSegmentHeader5(daf::DAF, desc::DAFSegmentDescriptor)

    iaa = initial_address(desc)
    faa = final_address(desc)

    i0 = 8*(faa-2)

    # Gravitational constant 
    GM = get_float(array(daf), i0, endian(daf))

    # Number of states stored in the record
    n = Int(get_float(array(daf), i0 + 8, endian(daf)))

    # Number of epoch directories 
    ndirs = n ÷ 100

    # Initial address for the epoch table (after all the state records)
    etid = 8*(iaa - 1) + 48*n # 8*6*n

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

    # Pre-compute the pairs required for the Stumpff functions series evaluation.
    # The truncation degree is set to 11, similarly to SPICE.
    deg = 11  
    
    p = zeros(2(deg-1))
    @inbounds for j in eachindex(p)
        p[j] = 1/(j*(j+1))
    end

    SPKSegmentHeader5(GM, n, ndirs, etid, epochs, iaa, p)

end

""" 
    SPKSegmentCache5(head::SPKSegmentHeader5)

Initialise the cache for an SPK segment of type 5.
"""
function SPKSegmentCache5(head::SPKSegmentHeader5) 
    SPKSegmentCache5(
        TwoBodyUniversalCache(head.GM), 
        TwoBodyUniversalCache(head.GM),
        MVector{2, Float64}(zeros(2)),
        -2
    )
end

""" 
    SPKSegmentType2(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 5.
"""
function SPKSegmentType5(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader5(daf, desc)
    caches = [SPKSegmentCache5(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType5(header, caches)

end

@inline spk_field(::SPKSegmentType5) = SPK_SEGMENTLIST_MAPPING[5]

function spk_vector3(daf::DAF, seg::SPKSegmentType5, time::Number) 

    head = header(seg)
    data = cache(seg)

    # Retrieve the state coefficients
    index = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, index)

    # Retrieve the bracketing epochs
    @inbounds t1, t2 = data.epochs[1], data.epochs[2] 

    # Compute the propagation intervals
    Δt1 = time - t1
    Δt2 = time - t2

    # Propagate the two states
    p1, _ = propagate_twobody(data.c1, Δt1, head.pairs)
    p2, _ = propagate_twobody(data.c2, Δt2, head.pairs)

    # Compute the weighting coefficient
    k = π/(t2 - t1)
    t = k*(time - t1)

    w = 1/2*(1 + cos(t))

    # Interpolate the position
    return w*p1 + (1-w)*p2

end

function spk_vector6(daf::DAF, seg::SPKSegmentType5, time::Number)

    head = header(seg)
    data = cache(seg)
    
    # Retrieve the state coefficients
    index = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, index)

    # Retrieve the bracketing epochs
    @inbounds t1, t2 = data.epochs[1], data.epochs[2] 

    # Compute the propagation intervals
    Δt1 = time - t1
    Δt2 = time - t2

    # Propagate the two states
    p1, v1 = propagate_twobody(data.c1, Δt1, head.pairs)
    p2, v2 = propagate_twobody(data.c2, Δt2, head.pairs)

    # Compute the weighting coefficient and its derivative 
    k = π/(t2 - t1)
    t = k*(time - t1)

    w = 1/2*(1 + cos(t))
    dw = -k/2*sin(t)

    # Interpolate the position and velocity
    pos = w*p1 + (1 - w)*p2
    vel = w*v1 + (1 - w)*v2 + dw*(p1 - p2) 

    return vcat(pos, vel)

end

function spk_vector9(::DAF, ::SPKSegmentType5, ::Number)
    throw(jEph.EphemerisError("Order ≥ 2 cannot be computed on segments of type 5."))
end

function spk_vector12(::DAF, ::SPKSegmentType5, ::Number)
    throw(jEph.EphemerisError("Order ≥ 2 on segments of type 5."))
end


"""
    find_logical_record(daf::DAF, head::SPKSegmentHeader5, time::Number)
"""
function find_logical_record(daf::DAF, head::SPKSegmentHeader5, time::Number)

    # We still need to differentiate the two cases (with or without epoch directory), 
    # because in the first we are reading the final epochs from the header, in the second 
    # case, we are gradually loading them as the header contains the directory epochs! 

    # The index of the first epoch >= than time is computed in a 0-index notation 
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
    # or equal to the requested time.

    return index

end

"""
    get_coefficients!(daf::DAF, head::SPKSegmentHeader5, cache::SPKSegmentCache5, index::Int)
"""
function get_coefficients!(daf::DAF, head::SPKSegmentHeader5, cache::SPKSegmentCache5, 
    index::Int)

    # Pre-process the index to always load two states, if index is 0, it is moved to 1, 
    first = max(0, index-1)

    # Check whether the coefficients for this record are already loaded
    first == cache.id && return nothing

    # Address of desired logical record 
    i0 = 8*(head.iaa - 1) + 48*first # 48 = 6 (states) * 8 bytes

    # Update both caches
    @inbounds for j = 1:3 
        cache.c1.pos[j] = get_float(array(daf), i0 + 8*(j-1), endian(daf)) 
        cache.c1.vel[j] = get_float(array(daf), i0 + 8*(j+2), endian(daf))

        cache.c2.pos[j] = get_float(array(daf), i0 + 48 + 8*(j-1), endian(daf)) 
        cache.c2.vel[j] = get_float(array(daf), i0 + 48 + 8*(j+2), endian(daf))
    end

    update_cache!(cache.c1)
    update_cache!(cache.c2)

    cache.id = first 

    # Here we load the epoch values 
    if head.ndirs == 0 
        # Extract the epochs directly from the ones stored in the header since there is 
        # no directory epoch
        @inbounds for j = 1:2
            cache.epochs[j] = head.epochs[first+j]
        end

    else 
        # Read the epochs from the binary file 
        i0 = head.etid + 8*first
        @inbounds for j = 1:2
            cache.epochs[j] = get_float(array(daf), i0 + 8*(j-1), endian(daf))
        end
    end

    nothing 
end
