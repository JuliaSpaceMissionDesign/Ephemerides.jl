
struct SPKSegmentType21
    n::Int                      # number of records in segment 
    ndirs::Int
    maxdim::Int                 # Difference wrt to type 1 
    dlsize::Int                 # Difference wrt to type 1 = 4*maxdim + 11
    tl::Vector{Float64}
    g::Vector{Float64}
    refpos::Vector{Float64}
    refvel::Vector{Float64}
    dt::Matrix{Float64}
    kqmax::Vector{Int}
    kq::Vector{Int}
    fc::Vector{Float64}
    wc::Vector{Float64}
    w::Vector{Float64}
end

function SPKSegmentType21(array, head::SPKSegmentHeader, lend::Bool)

    # number of records in segment
    i0 = 8*(head.faa-2)
    
    maxdim = Int(get_float(array, i0, lend))
    dlsize = 4*maxdim + 11

    # get number of directories in this record
    n = Int(get_float(array, i0 + 8, lend))
    ndirs = n ÷ 100

    SPKSegmentType21(
        n, ndirs, maxdim, dlsize, Int[0.0], zeros(maxdim), zeros(3), zeros(3), 
        zeros(maxdim, 3), Int[0.0], zeros(Int, 3), zeros(14), zeros(13), zeros(17)
    )
end

function find_logical_record(array, head::SPKSegmentHeader, seg::SPKSegmentType21, 
            lend::Bool, time::Number)

    i0 = 8*(head.faa-1)

    # search target epoch in epoch directories
    idx, found = 0, false 
    while !found && idx <= seg.ndirs 
        idx += 1

        epoch_dir = get_float(array, i0 - 8*(seg.ndirs - idx + 1), lend)
        if epoch_dir >= time 
            found = true 
        end

    end

    # Compute epoch table start and stop indexes! 
    if found 
        start_idx = (idx - 1) * 100 + 1
        stop_idx = 100*idx
    else
        start_idx = 100*seg.ndirs + 1
        stop_idx = seg.n
    end

    # GET LOGICAL RECORD WITHIN THAT DIRECTORY

    # i0 per epoch table (dopo segmenti record)
    i0 = 8*(head.iaa - 1) + 8*seg.dlsize*seg.n

    # i1 = i0 + 8*(start_idx - 1)

    record_idx, found = start_idx-1, false
    while record_idx <= stop_idx && !found 
        record_idx += 1 
        epoch_time = get_float(array, i0 + 8*(record_idx-1), lend)
        if epoch_time >= time 
            found = true 
        end
    end

    return found ? record_idx : stop_idx 

    # upper_boundary = get_float(array, i0 + 8*(record_idx-1), header.lend)
    # lower_boundary = record_idx != 1 ? 
                    # get_float(array, i0 + 8*(record_idx-2), header.lend) : 
                    # spkhead.tstart


end

function get_coefficients!(array, head::SPKSegmentHeader, seg::SPKSegmentType21, 
            index::Integer, lend::Bool)

    i0 = 8*(head.iaa - 1) + 8*seg.dlsize*(index-1)
    seg.tl[1] = get_float(array, i0, lend)
    
    for k = 1:seg.maxdim 
        seg.g[k] = get_float(array, i0 + 8k, lend)    
    end
    
    i2 = i0 + 8*(seg.maxdim+1)
    for k = 1:3 
        seg.refpos[k] = get_float(array, i2 + 16*(k-1), lend)
        seg.refvel[k] = get_float(array, i2 + 16k - 8, lend)
    end 
    
    i3 = i2 + 48

    for k = 1:3 
        for p = 1:seg.maxdim
            seg.dt[p, k] = get_float(array, i3 + 8*(p-1) + 8*seg.maxdim*(k-1), lend)
        end
    end
    
    i4 = i3 + 3*seg.maxdim*8
    seg.kqmax[1] = Int(get_float(array, i4, lend))
    
    for j = 1:3 
        seg.kq[j] = get_float(array, i4 + 8j, lend)
    end

end

function compute_mda3(seg::SPKSegmentType21, time::Number)
    
    Δ = time - seg.tl[1]
    tp = Δ
    mq2 = seg.kqmax[1] - 2
    ks = seg.kqmax[1] - 1
    
    seg.fc[1] = 1.0
    for j = 1:mq2
        seg.fc[j+1] = tp / seg.g[j] 
        seg.wc[j] = Δ / seg.g[j]
        tp = Δ + seg.g[j]
    end
    
    # compute inverse coefficients 
    for j = 1:seg.kqmax[1]
        seg.w[j] = 1.0 / j
    end
    
    jx = 0
    ks1 = ks - 1
    
    while ks >= 2 
        jx += 1 
    
        for j = 1:jx 
            seg.w[j+ks] = seg.fc[j+1]*seg.w[j+ks1] - seg.wc[j] * seg.w[j+ks] 
        end
    
        ks = ks1 
        ks1 -= 1
    end

    # Compute x 
    psum = 0.0 
    @simd for j = seg.kq[1]:-1:1 
        psum += seg.dt[j, 1]*seg.w[j+ks]
    end

    x = seg.refpos[1] + Δ*(seg.refvel[1] + Δ*psum)

    # Compute y
    psum = 0.0 
    @simd for j = seg.kq[2]:-1:1 
        psum += seg.dt[j, 2]*seg.w[j+ks]
    end 

    y = seg.refpos[2] + Δ*(seg.refvel[2] + Δ*psum)

    # Compute z 
    psum = 0.0 
    @simd for j = seg.kq[3]:-1:1 
        psum += seg.dt[j, 3]*seg.w[j+ks]
    end 

    z = seg.refpos[3] + Δ*(seg.refvel[3] + Δ*psum)

    return SA[x, y, z]
end

function vector3(array, daf::DAFHeader, head::SPKSegmentHeader, 
    seg::SPKSegmentType21, time::Number)

    index = find_logical_record(array, head, seg, daf.lend, time)
    get_coefficients!(array, head, seg, index, daf.lend)

    # Compute position! 
    compute_mda3(seg, time)

end

