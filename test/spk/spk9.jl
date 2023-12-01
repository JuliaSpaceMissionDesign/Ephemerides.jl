
test_dir = artifact"testdata"
DJ2000 = 2451545

@testset "SPK Type 9" verbose=true begin 

    axes = Dict(1 => "J2000", 17 => "ECLIPJ2000")

    # The first kernel has a directory epoch, the second does not
    kernels = [joinpath(test_dir, "spk9_ex1.bsp"),  # even, epoch dir
               joinpath(test_dir, "spk9_ex2.bsp"),  # even, no epoch dir 
               joinpath(test_dir, "spk9_ex3.bsp"),  # uneven, epoch dir
               joinpath(test_dir, "spk9_ex4.bsp")]  # uneven no epoch dir 

    for kernel in kernels 

        ephj = EphemerisProvider(kernel);
        ephc = CalcephProvider(kernel);
        furnsh(kernel)
        
        desc = ephj.files[1].desc[1]
        t1j, t2j = desc.tstart, desc.tend

        head = ephj.files[1].seglist.spk9[1].head

        cid = Int(desc.cid)
        tid = Int(desc.tid)

        # # Test values 
        # yc1 = zeros(3);
        # yc2 = zeros(6);
        # yc3 = zeros(9);
        # yc4 = zeros(12);

        ep = t1j:1:t2j
        for j in 1:3000
            
            if iseven(j) 
                cid, tid = tid, cid 
            end

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
        
            # jEphem.ephem_compute!(yc1, ephc, DJ2000, tc, tid, cid, 0);
            # jEphem.ephem_compute!(yc2, ephc, DJ2000, tc, tid, cid, 1);
            # jEphem.ephem_compute!(yc3, ephc, DJ2000, tc, tid, cid, 2);
            # jEphem.ephem_compute!(yc4, ephc, DJ2000, tc, tid, cid, 3);

            ys1 = spkpos("$tid", tj, axes[desc.axesid], "NONE", "$cid")[1]
            ys2 = spkezr("$tid", tj, axes[desc.axesid], "NONE", "$cid")[1]
    
            # Test against CALCEPH (removed because erroneous)
            # @test yj1 ≈ yc1 atol=1e-9 rtol=1e-9
            # @test yj2 ≈ yc2 atol=1e-9 rtol=1e-9
            # @test yj3 ≈ yc3 atol=1e-9 rtol=1e-9
            # @test yj4 ≈ yc4 atol=1e-9 rtol=1e-9

            # Test against SPICE
            @test yj1 ≈ ys1 atol=1e-13 rtol=1e-14
            @test yj2 ≈ ys2 atol=1e-13 rtol=1e-14

            # Test if AUTODIFF works 
            # Position (doesn't work exactly for SPK type 9 because the coefficients are 
            # different, so we have reduced tolerances)
            # @test D¹(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[4:6] atol=1e-5 rtol=1e-5
            # @test D²(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[7:9] atol=1e-5 rtol=1e-5
            # @test D³(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[10:12] atol=1e-5 rtol=1e-5

            # Velocity (these ensure acceleration is computed correctly)
            @test D¹(t->ephem_vector6(ephj, cid, tid, t), tj)[4:6] ≈ yj4[7:9] atol=1e-9 rtol=1e-13
            @test D²(t->ephem_vector6(ephj, cid, tid, t), tj)[4:6] ≈ yj4[10:12] atol=1e-9 rtol=1e-13

            # Acceleration (there ensure the jerk is computed correctly)
            @test D¹(t->ephem_vector9(ephj, cid, tid, t), tj)[4:9] ≈ yj4[7:12] atol=1e-9 rtol=1e-13

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
