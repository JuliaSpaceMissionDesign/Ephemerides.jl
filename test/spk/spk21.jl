
test_dir = artifact"testdata"
DJ2000 = 2451545

@testset "SPK Type 21" verbose=true begin 
    
    # These kernels are tested against SPICE because CALCEPH performs erroneous 
    # computations on SPK types 21 (they differ from those of SPICE)

    # The first kernel has no epoch directories, whereas the second one does
    kernels = [joinpath(test_dir, "example1spk_seg21.bsp"), 
               joinpath(test_dir, "example2spk_seg21.bsp")]

    yc1 = zeros(3)

    for kernel in kernels

        ephj = EphemerisProvider(kernel);
        furnsh(kernel)

        desc = ephj.files[1].desc[1]
        head = ephj.files[1].seglist.spk1[1].head

        # Center and target bodies 
        cid = Int(desc.cid)
        tid = Int(desc.tid)

        t1j, t2j = desc.tstart, desc.tend

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

            # Test with Julia
            yj1 = ephem_vector3(ephj, cid, tid, tj);
            yj2 = ephem_vector6(ephj, cid, tid, tj);

            # Test with SPICE
            ys1 = spkpos("$tid", tj, "J2000", "NONE", "$cid")[1]
            ys2 = spkezr("$tid", tj, "J2000", "NONE", "$cid")[1]

            @test yj1 ≈ ys1 atol=1e-13 rtol=1e-14
            @test yj2 ≈ ys2 atol=1e-13 rtol=1e-14

            # Test if AUTODIFF works 
            @test D¹(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj2[4:end] atol=1e-9 rtol=1e-13

            # TODO: implement the acceleration and jerk. This functions below cannot be tested!
            # D²(t->ephem_vector3(ephj, cid, tid, t), tj)
            # D³(t->ephem_vector3(ephj, cid, tid, t), tj)

            # D¹(t->ephem_vector6(ephj, cid, tid, t), tj)
            # D²(t->ephem_vector6(ephj, cid, tid, t), tj)
            # D³(t->ephem_vector6(ephj, cid, tid, t), tj)

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

end
