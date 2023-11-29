
test_dir = artifact"testdata"
DJ2000 = 2451545

@testset "SPK Type 5" verbose=true begin 

    # The first kernel doesnt contains a directory epoch, the second one does
    kernels = [
        joinpath(test_dir, "example1spk_seg5.bsp"), 
        joinpath(test_dir, "example2spk_seg5.bsp")
    ]

    for kernel in kernels 
    
        ephj = EphemerisProvider(kernel);
        furnsh(kernel)

        desc = ephj.files[1].desc[1]

        t1j, t2j = desc.tstart, desc.tend

        head = ephj.files[1].seglist.spk5[1].head

        cid = Int(desc.cid)
        tid = Int(desc.tid)
        
        # Check errors 
        @test_throws jEphem.EphemerisError ephem_vector9(ephj, cid, tid, t1j)
        @test_throws jEphem.EphemerisError ephem_vector12(ephj, cid, tid, t1j)  

        # Test values 
        yc1 = zeros(3);
        yc2 = zeros(6);

        ep = t1j:0.1:t2j
        for j in 1:2000
            
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
                tj = min(t2j, max(t1j, rand(head.epochs)))
            elseif j < 500
                # Test directory handling close to the borders
                tj = min(t2j, max(rand(head.epochs) + randn(), t1j))
            else 
                tj = rand(ep)
            end
            
            tc = tj/86400

            yj1 = ephem_vector3(ephj, cid, tid, tj);
            yj2 = ephem_vector6(ephj, cid, tid, tj);

            # jEphem.ephem_compute!(yc1, ephc, DJ2000, tc, tid, cid, 0);
            # jEphem.ephem_compute!(yc2, ephc, DJ2000, tc, tid, cid, 1);
            
            ys1 = spkpos("$tid", tj, "J2000", "NONE", "$cid")[1]
            ys2 = spkezr("$tid", tj, "J2000", "NONE", "$cid")[1]

            # Test against CALCEPH
            # @test yj1 ≈ yc1 atol=1e-6 rtol=1e-8
            # @test yj2 ≈ yc2 atol=1e-6 rtol=1e-8

            # Test against SPICE
            @test yj1 ≈ ys1 atol=1e-14 rtol=1e-15
            @test yj2 ≈ ys2 atol=1e-14 rtol=1e-15

            # For these segments autodiff does not have any physical meaning!

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


