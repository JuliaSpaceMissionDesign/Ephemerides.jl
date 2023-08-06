
""" 
    SPKSegmentHeader1(daf::DAF, desc::DAFSegmentDescriptor)
"""
function SPKSegmentHeader1(daf::DAF, desc::DAFSegmentDescriptor)

    # Number of records in segment
    n = Int(get_float(daf.array, 8*(desc.faa-1), daf.header.lend))
    SPKSegmentHeader1(n, n ÷ 100)

end

""" 
    SPKSegmentCache1()
"""
function SPKSegmentCache1()
    SPKSegmentCache1(
        Int[0.0], zeros(15), zeros(3), zeros(3), zeros(15, 3), 
        Int[0.0], zeros(Int, 3), zeros(14), zeros(13), zeros(17)
    )
end

""" 
    SPKSegmentType1(daf::DAF, desc::DAFSegmentDescriptor)
"""
function SPKSegmentType1(daf::DAF, desc::DAFSegmentDescriptor)

    SPKSegmentType1(
        SPKSegmentHeader1(daf, desc), 
        SPKSegmentCache1()
    )

end

@inline spk_field(spk::SPKSegmentType1) = 1

