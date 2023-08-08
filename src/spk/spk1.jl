
""" 
    SPKSegmentHeader1(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 1.
"""
function SPKSegmentHeader1(daf::DAF, desc::DAFSegmentDescriptor)

    # Number of records in segment
    n = Int(get_float(array(daf), 8*(final_address(desc)-1), endian(daf)))
    SPKSegmentHeader1(n, n ÷ 100)

end

""" 
    SPKSegmentCache1()

Initialise the cache for an SPK segment of type 1.
"""
function SPKSegmentCache1()
    SPKSegmentCache1(
        Int[0.0], zeros(15), zeros(3), zeros(3), zeros(15, 3), 
        Int[0.0], zeros(Int, 3), zeros(14), zeros(13), zeros(17)
    )
end

""" 
    SPKSegmentType1(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 1.
"""
function SPKSegmentType1(daf::DAF, desc::DAFSegmentDescriptor)

    SPKSegmentType1(
        SPKSegmentHeader1(daf, desc), 
        SPKSegmentCache1()
    )

end

@inline spk_field(::SPKSegmentType1) = SPK_SEGMENTLIST_MAPPING[1]

function spk_vector3(daf::DAF, seg::SPKSegmentType1, desc::DAFSegmentDescriptor, time::Number) 
    # TODO:
end


function spk_vector6(daf::DAF, seg::SPKSegmentType1, desc::DAFSegmentDescriptor, time::Number)
    # TODO:
end

function spk_vector9(::DAF, ::SPKSegmentType1, ::DAFSegmentDescriptor, ::Number)
    throw(EphemerisError("Order ≥ 2 cannot be computed on segments of type 1 and 21."))
end

function spk_vector12(::DAF, ::SPKSegmentType1, ::DAFSegmentDescriptor, ::Number)
    throw(EphemerisError("Order ≥ 2 on segments of type 1 and 21."))
end
