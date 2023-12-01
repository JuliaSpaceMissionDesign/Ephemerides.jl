
test_dir = artifact"testdata"
DJ2000 = 2451545

@testset "SPK Type 17" verbose=true begin 

    # The first kernel contains an odd-sized window, the second one an even-size window
    kernel = joinpath(test_dir, "spk17_ex1.bsp")
    
    ephj = EphemerisProvider(kernel);
    # ephc = CalcephProvider(kernel);
    furnsh(kernel)

    for desc in ephj.files[1].desc

        t1j, t2j = desc.tstart, desc.tend

        cid = Int(desc.cid)
        tid = Int(desc.tid)
    
        # Check errors 
        @test_throws jEphem.EphemerisError ephem_vector9(ephj, cid, tid, t1j)
        @test_throws jEphem.EphemerisError ephem_vector12(ephj, cid, tid, t1j)  

        # Test values 
        yc1 = zeros(3);
        yc2 = zeros(6);

        ep = t1j:1:t2j
        for j in 1:2000
            
            if iseven(j) 
                cid, tid = tid, cid 
            end

            tj = j == 1 ? t1j : (j == 2 ? t2j : rand(ep))
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
            @test yj1 ≈ ys1 atol=1e-10 rtol=1e-10
            @test yj2 ≈ ys2 atol=1e-10 rtol=1e-10

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

    end 

    kclear()


end;


