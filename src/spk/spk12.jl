
""" 
    SPKSegmentType12(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 12.
"""
function SPKSegmentType12(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader8(daf, desc)
    caches = [SPKSegmentCache8(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType12(header, caches)

end

@inline spk_field(::SPKSegmentType12) = SPK_SEGMENTLIST_MAPPING[12]


function spk_vector3(daf::DAF, seg::SPKSegmentType12, time::Number) 

    index = find_logical_record(header(seg), time)
    # get_coefficients!(daf, header(seg), cache(seg), index)

    # TODO: implement me
    return SA[0, 0, 0]

end

function spk_vector6(daf::DAF, seg::SPKSegmentType12, time::Number)

    index = find_logical_record(header(seg), time)
    # get_coefficients!(daf, header(seg), cache(seg), index)

    # TODO: implement me 
    return SA[0, 0, 0, 0, 0, 0]

end

function spk_vector9(daf::DAF, seg::SPKSegmentType12, time::Number)

    index = find_logical_record(header(seg), time)
    # get_coefficients!(daf, header(seg), cache(seg), index)

    # TODO: implement me
    return SA[0, 0, 0, 0, 0, 0, 0, 0, 0]

end

function spk_vector12(daf::DAF, seg::SPKSegmentType12, time::Number)

    index = find_logical_record(header(seg), time)
    # get_coefficients!(daf, header(seg), cache(seg), index)

    # TODO: implement me 
    return SA[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
end