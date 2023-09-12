
test_dir = artifact"testdata"

DJ2000 = 2451545

@testset "SPK Type 18" verbose=true begin 

    # Each kernel contains segments both with and without epoch directories. 
    # Different descriptors are tested to include epoch directories and a number 
    # of points smaller than the interpolating window

    kernels = [joinpath(test_dir, "example1spk_seg18.bsp"),  # subtype = 0
               joinpath(test_dir, "example2spk_seg18.bsp"),  # subtype = 0
               joinpath(test_dir, "example3spk_seg18.bsp")]  # subtype = 0, with n < N

    # Desired descriptors 
    descidj = [(1, 2), (2, 3, 4), (1,)]

    for (ik, kernel) in enumerate(kernels)

        ephj = EphemerisProvider(kernel);
        furnsh(kernel)

        for idj in descidj[ik]
            desc = ephj.files[1].desc[idj]
            t1j, t2j = desc.tstart, desc.tend

            head = ephj.files[1].seglist.spk18[idj].head

            cid = Int(desc.cid)
            tid = Int(desc.tid)

            ep = t1j:1:t2j
            for j in 1:3000

                if j == 1 
                    # Test initial time
                    tj = t1j 
                elseif j == 2
                    # Test final time
                    tj = t2j 
                elseif j < 100
                    # Test values at the directory epochs
                    tj = rand(head.epochs)
                elseif j < 500
                    # Test directory handling close to the borders
                    tj = min(t2j, max(rand(head.epochs) + randn(), t1j))
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
                # Position (doesn't work exactly because pos and vel have different set of coefficients)
                # @test D¹(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[4:6] atol=1e-9 rtol=1e-13
                # @test D²(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[7:9] atol=1e-9 rtol=1e-13
                # @test D³(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[10:12] atol=1e-9 rtol=1e-13

                # Velocity (these test also the acceleration\jerk values that are computed)
                @test D¹(t->ephem_vector6(ephj, cid, tid, t), tj)[4:6] ≈ yj3[7:9] atol=1e-8 rtol=1e-9
                @test D²(t->ephem_vector6(ephj, cid, tid, t), tj)[4:6] ≈ yj4[10:12] atol=1e-8 rtol=1e-9

                # Acceleration  (these further test the jerk)
                @test D¹(t->ephem_vector9(ephj, cid, tid, t), tj)[4:9] ≈ yj4[7:12] atol=1e-8 rtol=1e-9

            end

        end

        kclear()

    end

end;


