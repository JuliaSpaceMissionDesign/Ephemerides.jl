
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
    order = Int(get_float(array(daf), i0 + 16, endian(daf)))
    N = order + 1

    # Number of states stores in the record 
    n = Int(get_float(array(daf), i0 + 24, endian(daf)))

    # Check even group size 
    iseven = N % 2 == 0

    SPKSegmentHeader8(tstart, tlen, order, N, n, iaa, iseven)
end

""" 
    SPKSegmentCache8(spkhead::SPKSegmentHeader2)

Initialise the cache for an SPK segment of type 8.
"""
function SPKSegmentCache8(header::SPKSegmentHeader8) 
    SPKSegmentCache8(
        zeros(header.N, 6), 
        zeros(header.N),
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

    # TODO: finish me
    return SA[time, time, time]

end


function spk_vector6(daf::DAF, seg::SPKSegmentType8, time::Number)

    # TODO: finish me
    return SA[time, time, time, 
              time, time, time]
end

function spk_vector9(daf::DAF, seg::SPKSegmentType8, time::Number)

    # TODO: finish me
    return SA[time, time, time, 
              time, time, time, 
              time, time, time]
end

function spk_vector12(daf::DAF, seg::SPKSegmentType8, time::Number)

    # TODO: finish me
    return SA[time, time, time, 
              time, time, time, 
              time, time, time, 
              time, time, time]
end


"""
    find_logical_record(head::SPKSegmentHeader8, time::Number)
"""
function find_logical_record(head::SPKSegmentHeader8, time::Number)

    Δt = time - seg.tstart 

    if head.iseven # even group size 
        low = Int(Δt ÷ head.tlen) + 1
        first = low - (head.N ÷ 2) + 1
    else # odd group size  
        near = round(Int, Δt/head.tlen) + 1
        first = near - head.order ÷ 2
    end

    index = min(max(1, first), head.n - head.order)
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
