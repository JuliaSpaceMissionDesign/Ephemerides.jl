
@testset "Properties" verbose=true begin 

    # Test segment boundaries
    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1.0, 1.6, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 1.6, 2, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 0.4, 1.0, 31008, 1, -1, 1, 1),
    ]
    
    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [0.4]
    @test te == [2]
    
    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1.0, 1.6, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 1.6, 2, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 0.4, 0.7, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 0.2, 2.8, 31008, 1, -1, 1, 1),
    ]
    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [0.2]
    @test te == [2.8]

    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1.0, 1.6, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 0.2, 0.8, 31008, 1, -1, 1, 1),
    ]

    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [0.2, 1.0]
    @test te == [0.8, 1.6]

end;