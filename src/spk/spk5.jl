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

    # Initial address for teh epoch table (after all the state records)
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
    
    SPKSegmentHeader5(GM, n, ndirs, etid, epochs)

end

""" 
    SPKSegmentCache5(head::SPKSegmentHeader2)

Initialise the cache for an SPK segment of type 5.
"""
function SPKSegmentCache5() 
    SPKSegmentCache5(
        MVector{6, Float64}(zeros(6)), 
        MVector{6, Float64}(zeros(6)),
        MVector{2, Float64}(zeros(2)),
        MVector(-1)
    )
end

""" 
    SPKSegmentType2(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 5.
"""
function SPKSegmentType5(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader5(daf, desc)
    caches = [SPKSegmentCache5() for _ in 1:Threads.nthreads()]

    SPKSegmentType5(header, caches)

end

@inline spk_field(::SPKSegmentType5) = SPK_SEGMENTLIST_MAPPING[5]

function spk_vector3(daf::DAF, seg::SPKSegmentType5, time::Number) 

    head = header(seg)
    data = cache(seg)

    x, y, z = zeros(3)
    return SA[x, y, z]

end


function spk_vector6(daf::DAF, seg::SPKSegmentType5, time::Number)

    head = header(seg)
    data = cache(seg)

    x, y, z, vx, vy, vz = zeros(6)
    return SA[x, y, z, vx, vy, vz]

end

function spk_vector9(::DAF, ::SPKSegmentType5, ::Number)
    throw(jEph.EphemerisError("Order ≥ 2 cannot be computed on segments of type 5."))
end

function spk_vector12(::DAF, ::SPKSegmentType5, ::Number)
    throw(jEph.EphemerisError("Order ≥ 2 on segments of type 5."))
end


