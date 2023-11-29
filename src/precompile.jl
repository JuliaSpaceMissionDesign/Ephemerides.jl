
PrecompileTools.@setup_workload begin

    test_dir = artifact"testdata"

    spks = [
        "example1spk_seg1.bsp",
        "example1spk_seg2.bsp",
        "example1spk_seg3.bsp",
        "example1spk_seg5.bsp",
        "example1spk_seg8.bsp",
        "example1spk_seg9.bsp",
        "example1spk_seg12.bsp",
        "example1spk_seg13.bsp",
        "example1spk_seg14.bsp",
        "example1spk_seg17.bsp",
        "example1spk_seg18.bsp",
        "example1spk_seg19.bsp",
        "example1spk_seg20.bsp",
        "example1spk_seg21.bsp",
    ]

    pcks = [
        "pa421.bpc",
    ]

    PrecompileTools.@compile_workload begin

        # Precompile SPK routines 
        for (j, file) in enumerate(spks) 
            kernel = joinpath(test_dir, file)
            ephj = EphemerisProvider(kernel)

            desc = ephj.files[1].desc[1]
            
            cid, tid = center(desc), target(desc)
            tj = desc.tend 

            _ = ephem_vector3(ephj, cid, tid, tj);
            _ = ephem_vector6(ephj, cid, tid, tj);

            if !(segment_type(desc) in (1, 5, 17, 21))
                _ = ephem_vector9(ephj, cid, tid, tj);
                _ = ephem_vector12(ephj, cid, tid, tj);
            end
        end

        # Precomile PCK routines 
        for file in pcks 
            kernel = joinpath(test_dir, file)
            ephj = EphemerisProvider(kernel)
    
            desc = ephj.files[1].desc[1]
    
            cid, tid = center(desc), target(desc)
            tj = desc.tend 
    
            _ = ephem_rotation3(ephj, cid, tid, tj);
            _ = ephem_rotation6(ephj, cid, tid, tj);
            _ = ephem_rotation9(ephj, cid, tid, tj);
            _ = ephem_rotation12(ephj, cid, tid, tj);
    
        end

    end

end