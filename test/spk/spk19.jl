
test_dir = artifact"testdata"

DJ2000 = 2451545

@testset "SPK Type 19" verbose=true begin 

    # Each kernel contains segments with and without epoch directories. 
    # The first kernel contains only subtype 2 mini-segments, whereas the second one 
    # has subtype 0 mini-segments. Both 0 and 1 are already well-tested in the spk18 
    # tests, so it is not strictly necessary to repeat them.

    # TODO: missing a test kernel with a boundary flag to select the first interval 
    
    kernels = [joinpath(test_dir, "example1spk_seg19.bsp"),  # subtype = 2 no directory
               joinpath(test_dir, "example2spk_seg19.bsp")]  # subtype = 0 has directory


    for (ik, kernel) in enumerate(kernels)

        ephj = EphemerisProvider(kernel);
        furnsh(kernel)

        desc = ephj.files[1].desc[1]
        t1j, t2j = desc.tstart, desc.tend

        head = ephj.files[1].seglist.spk19[1].head

        cid = Int(desc.cid)
        tid = Int(desc.tid)

        ep = t1j:1:t2j
        for j in 1:3000
            
            if iseven(j) 
                cid, tid = tid, cid 
            end

            if j == 1 
                # Test start time 
                tj = t1j 
            elseif j == 2 
                # Test end time
                tj = t2j 
            elseif j < 100
                # Test correct boundary flag handling
                tj = rand(head.times)
            elseif j < 500 
                # Test boundary flag handling close to the borders
                tj = min(t2j, max(rand(head.times) + randn(), t1j))
            else 
                tj = rand(ep)
            end

            tc = tj/86400

            yj1 = ephem_vector3(ephj, cid, tid, tj);
            yj2 = ephem_vector6(ephj, cid, tid, tj);
            yj3 = ephem_vector9(ephj, cid, tid, tj);
            yj4 = ephem_vector12(ephj, cid, tid, tj);
                    
            ys1 = spkpos("$tid", tj, "J2000", "NONE", "$cid")[1]
            ys2 = spkezr("$tid", tj, "J2000", "NONE", "$cid")[1]

            # Test against SPICE
            @test yj1 ≈ ys1 atol=1e-13 rtol=1e-14
            @test yj2 ≈ ys2 atol=1e-13 rtol=1e-14

            # Test if AUTODIFF works 
            # Position (doesn't work exactly because pos and vel might have different 
            # set of polynomial coefficients depending on the actual mini-segment subtype)
            # @test D¹(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[4:6] atol=1e-9 rtol=1e-13
            # @test D²(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[7:9] atol=1e-9 rtol=1e-13
            # @test D³(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[10:12] atol=1e-9 rtol=1e-13

            # When we are close to the segment boundaries, the accuracy reduces because 
            # of the lower order of the polynomials 

            # Velocity (these test also the acceleration\jerk values that are computed)
            @test D¹(t->ephem_vector6(ephj, cid, tid, t), tj)[4:6] ≈ yj3[7:9] atol=1e-5 rtol=1e-13
            @test D²(t->ephem_vector6(ephj, cid, tid, t), tj)[4:6] ≈ yj4[10:12] atol=1e-5 rtol=1e-13

            # Acceleration  (these further test the jerk)
            @test D¹(t->ephem_vector9(ephj, cid, tid, t), tj)[4:9] ≈ yj4[7:12] atol=1e-5 rtol=1e-13

        end

        # Thread-safe testing 
        tj = shuffle(collect(LinRange(t1j, t2j, 200)))

        pos = zeros(3, length(tj))
        for j in eachindex(tj)
            pos[:, j] =  ephem_vector3(ephj, cid, tid, tj[j])
        end
    
        pos_m = zeros(3, length(tj))
        Threads.@threads for j in eachindex(tj)
            pos_m[:, j] = ephem_vector3(ephj, cid, tid, tj[j])
        end

        @test pos ≈ pos_m atol=1e-14 rtol=1e-14

        kclear()

    end

end;


