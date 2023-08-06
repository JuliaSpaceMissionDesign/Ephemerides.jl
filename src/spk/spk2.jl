
""" 
    SPKSegmentHeader2(daf::DAF, desc::DAFSegmentDescriptor)
"""
function SPKSegmentHeader2(daf::DAF, desc::DAFSegmentDescriptor)

    i0 = 8*(desc.faa-4)

    tstart = get_float(daf.array, i0, daf.header.lend)
    tlen = get_float(daf.array, i0+8, daf.header.lend)

    # polynomial order 
    rsize = Int(get_float(daf.array, i0+16, daf.header.lend))

    # The order of the polynomial is actually = (order-1)
    order = (rsize - 2) ÷ 3

    # number of records
    n = Int(get_float(daf.array, i0+24, daf.header.lend))

    SPKSegmentHeader2(tstart, tlen, order, n)

end

""" 
    SPKSegmentCache2(spkhead::SPKSegmentHeader2)
"""
function SPKSegmentCache2(spkhead::SPKSegmentHeader2) 
    SPKSegmentCache2(zeros(3, spkhead.order), zeros(spkhead.order))
end

""" 
    SPKSegmentType2(daf::DAF, desc::DAFSegmentDescriptor)
"""
function SPKSegmentType2(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    head = SPKSegmentHeader2(daf, desc)
    cache = SPKSegmentCache2(head)

    SPKSegmentType2(head, cache)

end

@inline spk_field(spk::SPKSegmentType2) = 2


function ephem_vector3(daf::DAF, seg::SPKSegmentType2, desc::DAFSegmentDescriptor, time::Number) 

    index, tfrac = find_logical_record(seg, time)
    get_coefficients!(daf, seg, desc, index);
    chebyshev!(seg, tfrac)

    x, y, z = 0.0, 0.0, 0.0
    @inbounds @simd for i in eachindex(seg.cache.x)
        x += seg.cache.x[i]*seg.cache.A[1, i]
        y += seg.cache.x[i]*seg.cache.A[2, i]
        z += seg.cache.x[i]*seg.cache.A[3, i]
    end

    return SA[x, y, z]

end

function find_logical_record(seg::SPKSegmentType2, time::Number)

    idx, tfrac = divrem(time - seg.head.tstart, seg.head.tlen) 
    index = round(Int, idx)

    if index == seg.head.n 
        # This should only happen when time equals the final segment time
        index -= 1 
        tfrac = seg.head.tlen 
    end 

    return index, tfrac
end


function get_coefficients!(daf::DAF, seg::SPKSegmentType2, 
            desc::DAFSegmentDescriptor, index::Integer)

    ncomp = 3 # components (3 for pos, 6 for vel)

    # size of each logical record in bytes for spk type 2
    recsize = 8*(seg.head.order*ncomp+2)

    # address of desired logical record (skipping mid and radius because they are all = )
    k = 8*(desc.iaa-1) + recsize*index + 16

    @inbounds for j = 1:seg.head.order 
        for i = 1:ncomp
            seg.cache.A[i, j] = get_float(
                daf.array, k + 8*(j-1) + 8*(i-1)*seg.head.order, daf.header.lend
            )
        end
    end

end

function chebyshev!(seg::SPKSegmentType2, tfrac::Number)

    seg.cache.x[1] = 1.0 
    seg.cache.x[2] = 2*tfrac/seg.head.tlen - 1 
    @inbounds for i = 3:seg.head.order 
        seg.cache.x[i] = 2.0*seg.cache.x[2]*seg.cache.x[i-1] - seg.cache.x[i-2]
    end

end 