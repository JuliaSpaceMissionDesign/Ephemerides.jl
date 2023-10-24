
@testset "Properties" verbose=true begin 

    # Test segment boundaries
    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1.0, 1.6, 31008, 1, -1, 1, 1),
    ]

    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [1.0]
    @test te == [1.6]

    # Test segment concatenated at the start
    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1.0, 1.6, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 0.5, 1.0, 31008, 1, -1, 1, 1)
    ]
    
    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [0.5]
    @test te == [1.6]

    # Test segment before the first
    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1.0, 1.6, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 0.5, 0.9, 31008, 1, -1, 1, 1)
    ]
    
    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [0.5, 1.0]
    @test te == [0.9, 1.6]

    # Test segment starts before the first and ends in the middle (c != b)
    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1, 2., 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 3, 4, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 5, 6, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 0.5, 4.5, 31008, 1, -1, 1, 1),
    ]

    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [0.5, 5.0]
    @test te == [4.5, 6.0]


    # Test segment concatenated at the end
    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1.0, 1.6, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 1.6, 2, 31008, 1, -1, 1, 1)
    ]
    
    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [1.0]
    @test te == [2]

    # Test segment after the first
    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1.0, 1.6, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 1.7, 2, 31008, 1, -1, 1, 1)
    ]
    
    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [1.0, 1.7]
    @test te == [1.6, 2]

    # Test segment ends after the end and starts in the middle (c != a)
    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1, 2., 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 3, 4, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 5, 6, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 2.5, 6.5, 31008, 1, -1, 1, 1),
    ]

    ts, te = Ephemerides.get_segment_boundaries(desclist)

    @test ts == [1.0, 2.5]
    @test te == [2.0, 6.5]


    # random tests at the end and also inclusive!
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

    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1.0, 1.6, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 1.4, 1.5, 31008, 1, -1, 1, 1),
    ]

    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [1.0]
    @test te == [1.6]

    # Test last scenario with removal of 1 ts 
    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1.0, 1.6, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 2.0, 3.0, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 1.8, 2.5, 31008, 1, -1, 1, 1),
    ]

    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [1.0, 1.8]
    @test te == [1.6, 3.0]

    # Test last scenario with removal of 1 te 
    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1.0, 1.6, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 2.0, 3.0, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 1.4, 1.8, 31008, 1, -1, 1, 1),
    ]

    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [1.0, 2.0]
    @test te == [1.8, 3.0]

    # test with removal of multiple points
    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1.0, 1.5, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 1.9, 2.4, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 2.7, 3.0, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 1.7, 2.85, 31008, 1, -1, 1, 1),
    ]
    
    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [1.0, 1.7]
    @test te == [1.5, 3.0]

    # test with removal of multiple points
    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1.0, 1.5, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 1.9, 2.4, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 2.7, 3.0, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 1.4, 2.5, 31008, 1, -1, 1, 1),
    ]
    
    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [1.0, 2.7]
    @test te == [2.5, 3.0]

    # test with removal of multiple points
    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1.0, 1.5, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 1.9, 2.4, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 2.7, 3.0, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 1.4, 2.3, 31008, 1, -1, 1, 1),
    ]
    
    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [1.0, 2.7]
    @test te == [2.4, 3.0]
    
    desclist = [
        Ephemerides.DAFSegmentDescriptor(2, 1, 2., 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 3, 4, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 5, 6, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 7, 8, 31008, 1, -1, 1, 1),
        Ephemerides.DAFSegmentDescriptor(2, 2.5, 4.5, 31008, 1, -1, 1, 1),
    ]
    
    ts, te = Ephemerides.get_segment_boundaries(desclist)
    @test ts == [1.0, 2.5, 5.0, 7.0]
    @test te == [2.0, 4.5, 6.0, 8.0]
    
    # Test initial and final times on SPK/PCK records
    rec_pck = Ephemerides.EphemRecordPCK(10, 20, [1.0, 2.0], [1.6, 2.1])
    @test Ephemerides.initial_times(rec_pck) == [1.0, 2.0]
    @test Ephemerides.final_times(rec_pck) == [1.6, 2.1]

    rec_spk = Ephemerides.EphemRecordSPK(10, 20, 2, [1.0, 2.0], [1.6, 2.1])
    @test Ephemerides.initial_times(rec_spk) == [1.0, 2.0]
    @test Ephemerides.final_times(rec_spk) == [1.6, 2.1]

end;