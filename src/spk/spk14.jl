

# function SegmentMetaData(daf::DAF, i0::Int)

#     # Number of meta-data elements
#     nmeta = Int(get_float(array(daf), i0, endian(daf)))

#     if nmeta != 17 
#         throw(ErrorException("unsupported number of segment meta data items."))
#     end

#     # Offset of the packet data from the start of a packet DAF_RECORD_LENGTH
#     pktoff = Int(get_float(array(daf), i0 - 8, endian(daf)))

#     # Size of the packets 
#     pktsz = Int(get_float(array(daf), i0 - 16, endian(daf)))

#     # Number of items in the reserved area 
#     nrsv = Int(get_float(array(daf), i0 - 24, endian(daf)))

#     # Base address of the reserved area 
#     rsvbas = Int(get_float(array(daf), i0 - 32, endian(daf)))

#     # Number of packets 
#     npkt = Int(get_float(array(daf), i0 - 40, endian(daf)))

#     # Base address of the packets 
#     pktbas = Int(get_float(array(daf), i0 - 48, endian(daf)))

#     # Type of the packet directory 
#     pdrtyp = Int(get_float(array(daf), i0 - 56, endian(daf)))

#     # Number of items in the packet directory 
#     npdr = Int(get_float(array(daf), i0 - 64, endian(daf)))

#     # Base address of the packet directory 
#     pdrbas = Int(get_float(array(daf), i0 - 72, endian(daf)))

#     # Number of reference items 
#     nref = Int(get_float(array(daf), i0 - 80, endian(daf)))

#     # Base address of the reference items 
#     refbas = Int(get_float(array(daf), i0 - 88, endian(daf)))

#     # Type of the reference directory 
#     rdrtyp = Int(get_float(array(daf), i0 - 96, endian(daf)))

#     # Number of items in the reference directory 
#     nrdr = Int(get_float(array(daf), i0 - 104, endian(daf)))

#     # Base address of the reference directory 
#     rdrbas = Int(get_float(array(daf), i0 - 112, endian(daf)))

#     # Number of constants in the segment 
#     ncon = Int(get_float(array(daf), i0 - 120, endian(daf)))

#     # Base Address of the constants in a generic segment.
#     conbas = Int(get_float(array(daf), i0 - 128, endian(daf)))
    
# end


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
        throw(ErrorException("unsupported number of segment meta data items."))
    end

    # Size of the packets 
    pktsize = Int(get_float(array(daf), i0 - 16, endian(daf)))

    # Number of packets 
    n = Int(get_float(array(daf), i0 - 40, endian(daf)))

    # Number of directory epochs 
    ndirs = Int(get_float(array(daf), i0 - 104, endian(daf)))

    # Retrieve polynomial order 
    order = get_float(array(daf), 8*(iaa-1), endian(daf))

    # Number of coefficients of the polynomial 
    N = order + 1

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

    SPKSegmentHeader14(order, n, ndirs, epochs, etid, ptid, pktsize, N)
end

""" 
    SPKSegmentCache14(head::SPKSegmentHeader2)

Initialise the cache for an SPK segment of type 14.
"""
function SPKSegmentCache14(head::SPKSegmentHeader14) 
    SPKSegmentCache14(MVector(-1))
end

""" 
    SPKSegmentType14(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 14.
"""
function SPKSegmentType14(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader14(daf, desc)
    caches = [SPKSegmentCache14(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType14(header, caches)

end

@inline spk_field(::SPKSegmentType14) = SPK_SEGMENTLIST_MAPPING[14]


function spk_vector3(daf::DAF, seg::SPKSegmentType14, time::Number) 

    head = header(seg)
    data = cache(seg)

    # Retrieve coefficients 
    index = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, index)

    x, y, z = zeros(3)
    return SA[x, y, z]

end

function spk_vector6(daf::DAF, seg::SPKSegmentType14, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve coefficients 
    index = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, index)

    x, y, z, vx, vy, vz = zeros(6)
    return SA[x, y, z, vx, vy, vz]

end

function spk_vector9(daf::DAF, seg::SPKSegmentType14, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve coefficients 
    index = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, index)

    x, y, z, vx, vy, vz, ax, ay, az = zeros(9)
    return SA[x, y, z, vx, vy, vz, ax, ay, az]

end

function spk_vector12(daf::DAF, seg::SPKSegmentType14, time::Number)

    head = header(seg)
    data = cache(seg)

    # Retrieve coefficients 
    index = find_logical_record(daf, head, time)
    get_coefficients!(daf, head, data, index)

    x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz = zeros(12)
    return SA[x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz]
    
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

    # The output index is returned in a 0-index notation
    return index

end

"""
    get_coefficients!(daf::DAF, head::SPKSegmentHeader14, cache::SPKSegmentCache14, index::Int)
"""
function get_coefficients!(daf::DAF, head::SPKSegmentHeader14, cache::SPKSegmentCache14, 
            index::Int)

    # Check whether the coefficients for this record are already loaded
    index == cache.id[1] && return nothing
    cache.id[1] = index 

    # Address of desired packet record 
    k = head.ptid + 8*head.pktsize*index 

    # Retrieve the record mid point and radius 
    cache.p[1] = get_float(array(daf), k, endian(daf))
    cache.p[2] = get_float(array(daf), k + 8, endian(daf))
    cache.p[3] = 1/cache.p[2]

    k += 16 

    # # TODO: can we speed-up this part by casting the byte content into the array at once?
    @inbounds for j = 1:head.N 
        for i = 1:6 
            cache.A[i, j] = get_float(
                array(daf), k + 8*(j-1) + 8*(i-1)*head.N, endian(daf)
            )
        end
    end

    nothing

end