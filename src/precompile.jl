
PrecompileTools.@setup_workload begin

    test_dir = artifact"testdata"

    spks = [
        "spk1_ex1.bsp",
        "spk2_ex1.bsp",
        "spk3_ex1.bsp",
        "spk5_ex1.bsp",
        "spk8_ex1.bsp",
        "spk9_ex1.bsp",
        "spk12_ex1.bsp",
        "spk13_ex1.bsp",
        "spk14_ex1.bsp",
        "spk15_ex1.bsp",
        "spk17_ex1.bsp",
        "spk18_ex1.bsp",
        "spk19_ex1.bsp",
        "spk20_ex1.bsp",
        "spk21_ex1.bsp",
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

            if !(segment_type(desc) in (1, 5, 15, 17, 21))
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