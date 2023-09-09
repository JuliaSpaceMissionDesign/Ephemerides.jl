
""" 
    SPKSegmentHeader8(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 8.
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
    SPKSegmentCache8(spkhead::SPKSegmentHeader2)

Initialise the cache for an SPK segment of type 8.
"""
function SPKSegmentCache8(header::SPKSegmentHeader8) 

    if header.type == 8 
        wsize = header.N 
        d³size = 0 # Only for SPK 12 we need an array to store the the third derivatives
    else 
        wsize = 2*header.N
        d³size = wsize 
    end 

    SPKSegmentCache8(
        zeros(header.N, 6), 
        DiffCache(zeros(wsize)),
        DiffCache(zeros(wsize)),
        DiffCache(zeros(wsize)),
        DiffCache(zeros(d³size)),
        MVector(-1)
    )
end

""" 
    SPKSegmentType8(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 8.
"""
function SPKSegmentType8(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader8(daf, desc)
    caches = [SPKSegmentCache8(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType8(header, caches)

end

@inline spk_field(::SPKSegmentType8) = SPK_SEGMENTLIST_MAPPING[8]


function spk_vector3(daf::DAF, seg::SPKSegmentType8, time::Number) 

    index = find_logical_record(header(seg), time)
    get_coefficients!(daf, header(seg), cache(seg), index)

    tbegin = header(seg).tstart + (index-1)*header(seg).tlen 
    Δt = (time - tbegin)/header(seg).tlen + 1

    x = lagrange(cache(seg), Δt, 1, header(seg).N)
    y = lagrange(cache(seg), Δt, 2, header(seg).N)
    z = lagrange(cache(seg), Δt, 3, header(seg).N)

    return SA[x, y, z]

end

function spk_vector6(daf::DAF, seg::SPKSegmentType8, time::Number)

    index = find_logical_record(header(seg), time)
    get_coefficients!(daf, header(seg), cache(seg), index)

    tbegin = header(seg).tstart + (index-1)*header(seg).tlen 
    Δt = (time - tbegin)/header(seg).tlen + 1

    x = lagrange(cache(seg), Δt, 1, header(seg).N)
    y = lagrange(cache(seg), Δt, 2, header(seg).N)
    z = lagrange(cache(seg), Δt, 3, header(seg).N)

    vx = lagrange(cache(seg), Δt, 4, header(seg).N)
    vy = lagrange(cache(seg), Δt, 5, header(seg).N)
    vz = lagrange(cache(seg), Δt, 6, header(seg).N)

    return SA[x, y, z, vx, vy, vz]

end

function spk_vector9(daf::DAF, seg::SPKSegmentType8, time::Number)

    index = find_logical_record(header(seg), time)
    get_coefficients!(daf, header(seg), cache(seg), index)

    tbegin = header(seg).tstart + (index-1)*header(seg).tlen 
    Δt = (time - tbegin)/header(seg).tlen + 1

    x = lagrange(cache(seg), Δt, 1, header(seg).N)
    y = lagrange(cache(seg), Δt, 2, header(seg).N)
    z = lagrange(cache(seg), Δt, 3, header(seg).N)
    
    vx, ax = ∂lagrange(cache(seg), Δt, 4, header(seg).N, header(seg).tlen)
    vy, ay = ∂lagrange(cache(seg), Δt, 5, header(seg).N, header(seg).tlen)
    vz, az = ∂lagrange(cache(seg), Δt, 6, header(seg).N, header(seg).tlen)

    return SA[x, y, z, vx, vy, vz, ax, ay, az]

end

function spk_vector12(daf::DAF, seg::SPKSegmentType8, time::Number)

    index = find_logical_record(header(seg), time)
    get_coefficients!(daf, header(seg), cache(seg), index)

    tbegin = header(seg).tstart + (index-1)*header(seg).tlen 
    Δt = (time - tbegin)/header(seg).tlen + 1

    x = lagrange(cache(seg), Δt, 1, header(seg).N)
    y = lagrange(cache(seg), Δt, 2, header(seg).N)
    z = lagrange(cache(seg), Δt, 3, header(seg).N)
    
    vx, ax, jx = ∂²lagrange(cache(seg), Δt, 4, header(seg).N, header(seg).tlen)
    vy, ay, jy = ∂²lagrange(cache(seg), Δt, 5, header(seg).N, header(seg).tlen)
    vz, az, jz = ∂²lagrange(cache(seg), Δt, 6, header(seg).N, header(seg).tlen)

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
    get_coefficients!(daf::DAF, head, cache, index::Int)
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

# This function leverages Neville's algorithm to recursively evaluate the 
# Lagrange polynomial's at the desired time.
function lagrange(cache::SPKSegmentCache8, x::Number, istate::Integer, N::Int)

    # This function is valid only for equally-spaced polynomials!
    # x is a re-work of the ascissa that starts at 1
    # istate is the index of the desired state

    work = get_tmp(cache.work, x)

    @inbounds for i = 1:N 
        work[i] = cache.states[i, istate]
    end

    # compute the lagrange polynomials using a recursive relationship
    @inbounds for j = 1:N-1 
        for i = 1:N-j 
            c1 = i + j - x
            c2 = x - i 

            work[i] = (c1*work[i] + c2*work[i+1])/j
        end
    end

    return work[1]

end 


function ∂lagrange(cache::SPKSegmentCache8, x::Number, istate::Integer, N::Int, Δt::Number)

    # This function is valid only for equally-spaced polynomials!

    work = get_tmp(cache.work, x)
    dwork = get_tmp(cache.dwork, x)

    @inbounds for i = 1:N 
        work[i] = cache.states[i, istate]
        dwork[i] = 0.0
    end

    # Precompute abscissa derivatives
    dc2 = 1/Δt
    dc1 = -dc2

    # compute the lagrange polynomials using a recursive relationship
    @inbounds for j = 1:N-1 
        for i = 1:N-j 
            c1 = i + j - x
            c2 = x - i 

            dwork[i] = (dc1*work[i]   + c1*dwork[i] + 
                        dc2*work[i+1] + c2*dwork[i+1])/j

            work[i] = (c1*work[i] + c2*work[i+1])/j

        end
    end

    return work[1], dwork[1]

end 

function ∂²lagrange(cache::SPKSegmentCache8, x::Number, istate::Integer, N::Int, Δt::Number)

    # This function is valid only for equally-spaced polynomials!
    
    work = get_tmp(cache.work, x)
    dwork = get_tmp(cache.dwork, x)
    ddwork = get_tmp(cache.ddwork, x)

    @inbounds for i = 1:N 
        work[i] = cache.states[i, istate]
        dwork[i] = 0.0
        ddwork[i] = 0.0
    end

    # Precompute abscissa derivatives
    dc2 = 1/Δt
    dc1 = -dc2

    # compute the lagrange polynomials using a recursive relationship
    @inbounds for j = 1:N-1 
        for i = 1:N-j 
            c1 = i + j - x
            c2 = x - i 

            ddwork[i] = (2*dc1*dwork[i]   + c1*ddwork[i] + 
                         2*dc2*dwork[i+1] + c2*ddwork[i+1])/j

            dwork[i] = (dc1*work[i]   + c1*dwork[i] + 
                        dc2*work[i+1] + c2*dwork[i+1])/j

            work[i] = (c1*work[i] + c2*work[i+1])/j

        end
    end

    return work[1], dwork[1], ddwork[1]

end 
