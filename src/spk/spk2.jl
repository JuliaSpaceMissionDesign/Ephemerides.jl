
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

