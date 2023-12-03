
""" 
    SPKSegmentHeader1(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 1.
"""
function SPKSegmentHeader1(daf::DAF, desc::DAFSegmentDescriptor)

    # Initial/final segment address
    iaa = initial_address(desc)
    faa = final_address(desc)

    # Number of records in segment
    n = Int(get_float(array(daf), 8*(faa-1), endian(daf)))

    # Number of directory epochs 
    ndir = n ÷ 100

    # Maximum dimension depending on segment type
    if segment_type(desc) == 1 
        maxdim = 15 
    elseif segment_type(desc) == 21
        maxdim = Int(get_float(array(daf), 8*(faa-2), endian(daf)))
    else 
        throw(jEph.EphemerisError("Invalid segment type."))
    end
    
    # Double precision numbers stored in each MDA record 
    recsize = max(71, 4*maxdim + 11)

    # Initial address for the epoch table (after all the segment records)
    etid = 8*(iaa - 1) + 8*recsize*n

    if ndir == 0
        # Load actual epochs
        epochs = zeros(n)
        @inbounds for j = 1:n 
            epochs[j] = get_float(array(daf), etid + 8*(j-1), endian(daf))
        end

    else 
        # Load directory epochs
        epochs = zeros(ndir)
        @inbounds for j = 1:ndir
            epochs[j] = get_float(array(daf), etid + 8*(n+j-1), endian(daf)) 
        end

    end

    SPKSegmentHeader1(n, ndir, epochs, iaa, etid, recsize, maxdim)

end

""" 
    SPKSegmentCache1(head::SPKSegmentHeader1)

Initialise the cache for an SPK segment of type 1.
"""
function SPKSegmentCache1(head::SPKSegmentHeader1)
    SPKSegmentCache1(
        0.0, zeros(head.maxdim), zeros(3), zeros(3), zeros(head.maxdim, 3), 
        0, zeros(Int, 3), -1, 
        DiffCache(zeros(head.maxdim-1)), DiffCache(zeros(head.maxdim-2)), 
        DiffCache(zeros(head.maxdim+2)), DiffCache(zeros(3))
    )
end

""" 
    SPKSegmentType1(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 1.
"""
function SPKSegmentType1(daf::DAF, desc::DAFSegmentDescriptor)

    header = SPKSegmentHeader1(daf, desc)
    caches = [SPKSegmentCache1(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType1(header, caches)

end

@inline spk_field(::SPKSegmentType1) = SPK_SEGMENTLIST_MAPPING[1]

function spk_vector3(daf::DAF, seg::SPKSegmentType1, time::Number) 
    
    # Find the logical record containing the MDA coefficients at `time`
    index = find_logical_record(daf, header(seg), time)
    
    # Retrieve the MDA coefficients 
    get_coefficients!(daf, header(seg), cache(seg), index)

    Δ = time - cache(seg).tl

    # Compute MDA position coefficients
    compute_mda_pos_coefficients!(cache(seg), Δ)
    return compute_mda_position(cache(seg), Δ)
    
end


function spk_vector6(daf::DAF, seg::SPKSegmentType1, time::Number)
    # Find the logical record containing the MDA coefficients at `time`
    index = find_logical_record(daf, header(seg), time)
    
    # Retrieve the MDA coefficients 
    get_coefficients!(daf, header(seg), cache(seg), index)
  
    Δ = time - cache(seg).tl

    # Compute MDA position coefficients
    compute_mda_pos_coefficients!(cache(seg), Δ)
    pos = compute_mda_position(cache(seg), Δ)

    # Compute MDA velocity coefficients
    compute_mda_vel_coefficients!(cache(seg), Δ)
    vel = compute_mda_velocity(cache(seg), Δ)

    return @inbounds SVector{6}(
        pos[1], pos[2], pos[3], 
        vel[1], vel[2], vel[3]
    )
end

function spk_vector9(::DAF, ::SPKSegmentType1, ::Number)
    throw(jEph.EphemerisError("Order ≥ 2 cannot be computed on segments of type 1 and 21."))
end

function spk_vector12(::DAF, ::SPKSegmentType1, ::Number)
    throw(jEph.EphemerisError("Order ≥ 2 on segments of type 1 and 21."))
end



"""
    find_logical_record(daf::DAF, head::SPKSegmentHeader1, time::Number)
"""
function find_logical_record(daf::DAF, head::SPKSegmentHeader1, time::Number)

    index = 0
    if head.ndirs == 0
        # Search through epoch table 
        while (index < head.n - 1 && head.epochs[index+1] < time)
            index += 1 
        end

        return index 
    end

    # First search through epoch directories
    subdir = 0
    while (subdir < head.ndirs && head.epochs[subdir+1] < time)
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

    return index 

end

"""
    get_coefficients!(daf::DAF, head::SPKSegmentHeader1, cache::SPKSegmentCache1, index::Int)
"""
function get_coefficients!(daf::DAF, head::SPKSegmentHeader1, cache::SPKSegmentCache1, 
            index::Int)
    
    # Check whether the coefficients for this record are already loaded
    index == cache.id && return nothing
    cache.id = index 

    # Initial record address 
    i0 = 8*(head.iaa - 1) + 8*head.recsize*index
    cache.tl = get_float(array(daf), i0, endian(daf))

    @inbounds for k = 1:head.maxdim 
        cache.g[k] = get_float(array(daf), i0 + 8k, endian(daf))    
    end
    
    i2 = i0 + 8*(head.maxdim+1)
    @inbounds  for k = 1:3 
        cache.refpos[k] = get_float(array(daf), i2 + 16*(k-1), endian(daf))
        cache.refvel[k] = get_float(array(daf), i2 + 16k - 8, endian(daf))
    end 
    
    i3 = i2 + 48

    @inbounds for k = 1:3 
        for p = 1:head.maxdim 
            cache.dt[p, k] = get_float(
                array(daf), i3 + 8*(p-1) + 8*head.maxdim*(k-1), endian(daf)
            )
        end
    end
    
    i4 = i3 + 3*head.maxdim*8
    cache.kqmax = Int(get_float(array(daf), i4, endian(daf)))
    
    @inbounds for j = 1:3 
        cache.kq[j] = get_float(array(daf), i4 + 8j, endian(daf))
    end

end

"""
    compute_mda_position(cache::SPKSegmentCache1, Δ::Number)
"""
function compute_mda_position(cache::SPKSegmentCache1, Δ::Number)
    
    @inbounds for k = 1:3 
        get_tmp(cache.vct, Δ)[k] = 0.0 
        @simd for j = cache.kq[k]:-1:1 
            get_tmp(cache.vct, Δ)[k] += cache.dt[j, k]*get_tmp(cache.w, Δ)[j+1]
        end
    end

    # At this point by definition of the above while loop, ks will always be = 1
    return SVector{3}(
        cache.refpos[1] + Δ*(cache.refvel[1] + Δ*get_tmp(cache.vct, Δ)[1]),
        cache.refpos[2] + Δ*(cache.refvel[2] + Δ*get_tmp(cache.vct, Δ)[2]),
        cache.refpos[3] + Δ*(cache.refvel[3] + Δ*get_tmp(cache.vct, Δ)[3])
    )

end

"""
    compute_mda_velocity(cache::SPKSegmentCache1, Δ::Number)
"""
function compute_mda_velocity(cache::SPKSegmentCache1, Δ::Number)

    @inbounds for k = 1:3 
        get_tmp(cache.vct, Δ)[k] = 0.0 
        @simd for j = cache.kq[k]:-1:1 
            get_tmp(cache.vct, Δ)[k] += cache.dt[j, k]*get_tmp(cache.w, Δ)[j]
        end
    end

    # At this point by definition of the above while loop, ks will always be = 0
    return SVector{3}(
        cache.refvel[1] + Δ*get_tmp(cache.vct, Δ)[1],
        cache.refvel[2] + Δ*get_tmp(cache.vct, Δ)[2],
        cache.refvel[3] + Δ*get_tmp(cache.vct, Δ)[3]
    )

end

# This function pre-computes the cache coefficients w. At the end of it, we will 
# always have ks = 1 and ks1 = 0 due to the definition of the last while loop
@inbounds function compute_mda_pos_coefficients!(cache::SPKSegmentCache1, Δ::Number)

    tp = Δ
    mq2 = cache.kqmax - 2
    ks = cache.kqmax - 1
    
    get_tmp(cache.fc, Δ)[1] = 1.0
    for j = 1:mq2
        get_tmp(cache.fc, Δ)[j+1] = tp / cache.g[j] 
        get_tmp(cache.wc, Δ)[j] = Δ / cache.g[j]
        tp = Δ + cache.g[j]
    end
    
    # compute inverse coefficients 
    for j = 1:cache.kqmax
        get_tmp(cache.w, Δ)[j] = 1.0 / j
    end
    
    jx = 0
    ks1 = ks - 1
    
    while ks >= 2 
        jx += 1 
    
        for j = 1:jx 
            get_tmp(cache.w, Δ)[j+ks] = get_tmp(cache.fc, Δ)[j+1]*get_tmp(cache.w, Δ)[j+ks1] - get_tmp(cache.wc, Δ)[j] * get_tmp(cache.w, Δ)[j+ks] 
        end
    
        ks = ks1 
        ks1 -= 1
    end

    nothing
end

# This function is such that the output ks = 0
@inbounds function compute_mda_vel_coefficients!(cache::SPKSegmentCache1, Δ::Number)

    for j = 1:cache.kqmax
        get_tmp(cache.w, Δ)[j+1] = get_tmp(cache.fc, Δ)[j+1]*get_tmp(cache.w, Δ)[j] - get_tmp(cache.wc, Δ)[j]*get_tmp(cache.w, Δ)[j+1]
    end
    
    nothing
end

