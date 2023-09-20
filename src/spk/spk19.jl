
""" 
    SPKSegmentHeader19(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 18 and 19.
"""
function SPKSegmentHeader19(daf::DAF, desc::DAFSegmentDescriptor)

    iaa = initial_address(desc)
    faa = final_address(desc)

    type = segment_type(desc)

    # SPK type 18 are handled as special cases of type 19 with a single mini-segment   
    # We're just trading away some small amount of memory because the constant type 18 
    # headers are duplicated inside each of the type 19 caches.
    type == 18 && return SPKSegmentHeader19(1, 0, Float64[], iaa, faa, 0, false, type)

    # If we have arrived here, the segment is a generic SPK type 19
    i0 = 8*(faa - 3)

    # Retrieve mini-segment stop pointer
    stop_point = Int(get_float(array(daf), i0, endian(daf)))

    # DAF byte address of the start interval times table
    etid = 8*(iaa - 1) + 8*(stop_point - 1)

    # Retrieve boundary flag
    flag = Int(get_float(array(daf), i0 + 8, endian(daf))) == 0

    # Retrieve number of mini-segments
    n = Int(get_float(array(daf), i0 + 16, endian(daf)))

    # DAF byte address of the start mini-segment pointers table
    ptid = i0 - 8n

    # Number of interval directories (if there are 100 intervals, no directory is stored)
    ndirs = (n-1) ÷ 100 

    if ndirs == 0 
        # Load actual interval start times
        times = zeros(n)
        @inbounds for j = 1:n
            times[j] = get_float(array(daf), etid + 8*(j-1), endian(daf))
        end
    else 
        times = zeros(ndirs)
        @inbounds for j = 1:ndirs
            times[j] = get_float(array(daf), etid + 8*(n+j), endian(daf))
        end
    end

    SPKSegmentHeader19(n, ndirs, times, iaa, etid, ptid, flag, type)
end

""" 
    SPKSegmentCache19(daf::DAF, head::SPKSegmentHeader2)

Initialise the cache for an SPK segment of type 19.
"""
function SPKSegmentCache19(daf::DAF, head::SPKSegmentHeader19) 

    # Need to look through all the mini-segments to retrieve the highest values 
    # of n and N to accomodate all possible vector and buffer sizes 

    n, N = -1, -1
    nbuff, buffsize, packetsize = -1, -1, -1

    for j = 1:head.n
        if j == head.n 
            # leverage the interval table byte address
            i0 = (head.type == 19 ? head.etid : 8*head.etid) - 8
        else
            # use the start pointer of the next mini-segment 
            pnt = Int(get_float(array(daf), head.ptid + 8*j, endian(daf)))
            i0 = 8*(head.iaa - 1) + 8*(pnt-2)
        end

        # Retrieve mini-segment subtype, window size and number of packets
        subtype = Int(get_float(array(daf), i0 - 16, endian(daf)))
        Nj = Int(get_float(array(daf), i0 - 8, endian(daf)))
        nj = Int(get_float(array(daf), i0, endian(daf)))

        # Update maximum window size
        N = max(N, Nj)

        # Update maximum epoch vector size
        ndirsj = (nj - 1) ÷ 100 
        n = max(n, ndirsj == 0 ? nj : ndirsj)    

        # Update the maximum packet size
        packetsize = max(packetsize, subtype == 0 ? 12 : 6)

        # Update the maximum number of buffers
        nbuff = max(nbuff, subtype == 2 ? 4 : 3)

        # Update the maximum buffer size 
        buffsize = max(buffsize, subtype == 1 ? Nj : 2*Nj)

    end

    # Create an empty SPK segment 18 mini-header and mini-cache using the largest buffer 
    minihead = SPKSegmentHeader18(n)
    minidata = SPKSegmentCache18(
        MVector(-1, -1, -1), 
        zeros(N), 
        zeros(N, packetsize), 
        InterpCache{Float64}(nbuff, buffsize)
    )

    id = MVector(-1)
    
    if head.type == 18 
        # Update the header with the single type 18 information 
        update_header!(minihead, daf, head.iaa, head.etid, head.type)
        @inbounds id[1] = 0 # Update the cached index
    end

    return SPKSegmentCache19(minihead, minidata, id)

end

""" 
    SPKSegmentType19(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 19.
"""
function SPKSegmentType19(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader19(daf, desc)
    caches = [SPKSegmentCache19(daf, header) for _ in 1:Threads.nthreads()]

    SPKSegmentType19(header, caches)

end

@inline spk_field(::SPKSegmentType19) = SPK_SEGMENTLIST_MAPPING[19]


function spk_vector3(daf::DAF, seg::SPKSegmentType19, time::Number) 

    head = header(seg)
    data = cache(seg)

    # Retrieve the index of the desired mini-segment
    index = find_minisegment(daf, head, time)

    # Load the mini-segment data 
    load_minisegment!(daf, head, data, index)

    # Interpolate inside the minisegment 
    return spk_vector3(daf, data.minihead, data.minidata, time)

end

function spk_vector6(daf::DAF, seg::SPKSegmentType19, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve the index of the desired mini-segment
    index = find_minisegment(daf, head, time)

    # Load the mini-segment data 
    load_minisegment!(daf, head, data, index)

    # Interpolate inside the minisegment 
    return spk_vector6(daf, data.minihead, data.minidata, time)

end

function spk_vector9(daf::DAF, seg::SPKSegmentType19, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve the index of the desired mini-segment
    index = find_minisegment(daf, head, time)

    # Load the mini-segment data 
    load_minisegment!(daf, head, data, index)

    # Interpolate inside the minisegment 
    return spk_vector9(daf, data.minihead, data.minidata, time)

end

function spk_vector12(daf::DAF, seg::SPKSegmentType19, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve the index of the desired mini-segment
    index = find_minisegment(daf, head, time)

    # Load the mini-segment data 
    load_minisegment!(daf, head, data, index)

    # Interpolate inside the minisegment 
    return spk_vector12(daf, data.minihead, data.minidata, time)

end


"""
    find_minirecord(daf::DAF, head::SPKSegmentHeader19, time::Number)
"""
function find_minisegment(daf::DAF, head::SPKSegmentHeader19, time::Number)

    index = 0
    if head.ndirs == 0
        # Search through interval table depending on the selected flag  
        @inbounds while (index < head.n - 1 && 
                        ((head.usefirst && head.times[index+2] <  time) || 
                        (!head.usefirst && head.times[index+2] <= time)))
            index += 1 
        end

    else 
        # Here, we first search through the interval directories, 
        # dependin on the boundary flag 
        subdir = 0
        @inbounds while (subdir < head.ndirs && 
                        ((head.usefirst && head.times[subdir+1] <  time) || 
                        (!head.usefirst && head.times[subdir+1] <= time)))
            subdir += 1 
        end

        index = subdir*100 - 1
        stop_idx = min(index + 100, head.n)

        # Now find the actual interval within the identified directory.
        found = false 
        while (index < stop_idx - 1 && !found)
            epoch = get_float(array(daf), head.etid + 8*(index+1), endian(daf))      

            if (epoch < time && head.usefirst) || (epoch <= time && !head.usefirst)
                index += 1 
            else 
                found = true 
            end 
        end

    end

    # At this point of the code we have found the index of the minisegment that we desire. 
    # Note that the index is in 0-based notation (i.e., the 1st miniseg has index 0)
    return index

end

"""
    load_minisegment!(daf::DAF, head::SPKSegmentHeader19, cache::SPKSegmentCache19, index::Int)
"""
function load_minisegment!(daf::DAF, head::SPKSegmentHeader19, cache::SPKSegmentCache19, 
            index::Int)

    # Check whether the desired mini-segment has already been loaded (i.e., equal to the 
    # one used in the previous call)
    @inbounds cache.id[1] == index && return nothing 
    @inbounds cache.id[1] = index

    # Retrive the initial and final addresses (1-based and not in byte) using the pointers
    iaa = head.iaa + Int(get_float(array(daf), head.ptid + 8*index, endian(daf))) - 1
    faa = head.iaa + Int(get_float(array(daf), head.ptid + 8*index + 8, endian(daf))) - 2

    # Update the mini-segment header with the new data 
    update_header!(cache.minihead, daf, iaa, faa, 19)
    
    # Reset the cache because the mini-segment has changed 
    reset_indexes!(cache.minidata)

    nothing

end 
