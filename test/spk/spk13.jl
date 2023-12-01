
test_dir = artifact"testdata"
DJ2000 = 2451545

@testset "SPK Type 13" verbose=true begin 

    # Each kernel contains segments both with and without epoch directories. 
    # Different descriptors are tested to include epoch directories and even/odd windows
    kernels = [joinpath(test_dir, "spk13_ex1.bsp"), 
               joinpath(test_dir, "spk13_ex2.bsp")]
    
    # Desired descriptors 
    descidj = [(1, 2), (3, 4)]

    for (ik, kernel) in enumerate(kernels)

        ephj = EphemerisProvider(kernel);
        ephc = CalcephProvider(kernel);
        furnsh(kernel)

        for (cdj, idj) in enumerate(descidj[ik])
            desc = ephj.files[1].desc[idj]
            head = ephj.files[1].seglist.spk9[idj].head

            t1j, t2j = desc.tstart, desc.tend

            cid = Int(desc.cid)
            tid = Int(desc.tid)

            # Test values 
            yc1 = zeros(3);
            yc2 = zeros(6);
            yc3 = zeros(9);
            yc4 = zeros(12);

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
            
                jEphem.ephem_compute!(yc1, ephc, DJ2000, tc, tid, cid, 0);
                jEphem.ephem_compute!(yc2, ephc, DJ2000, tc, tid, cid, 1);
                jEphem.ephem_compute!(yc3, ephc, DJ2000, tc, tid, cid, 2);
                jEphem.ephem_compute!(yc4, ephc, DJ2000, tc, tid, cid, 3);
                        
                ys1 = spkpos("$tid", tj, "J2000", "NONE", "$cid")[1]
                ys2 = spkezr("$tid", tj, "J2000", "NONE", "$cid")[1]

                # Test against CALCEPH
                @test yj1 ≈ yc1 atol=1e-9 rtol=1e-9
                @test yj2 ≈ yc2 atol=1e-9 rtol=1e-9
                @test yj3 ≈ yc3 atol=1e-9 rtol=1e-9
                @test yj4 ≈ yc4 atol=1e-9 rtol=1e-9

                # Test against SPICE
                @test yj1 ≈ ys1 atol=1e-13 rtol=1e-14
                @test yj2 ≈ ys2 atol=1e-13 rtol=1e-14

                # The accuracy of the interpolation derivatives reduces when the requested
                # time is close to the segment boundaries 
                if ik == 2 && cdj == 2
                    jatol, jrtol = 1e-5, 1e-9 
                else 
                    jatol, jrtol = 1e-5, 1e-13
                end

                # Test if AUTODIFF works 
                # Position
                @test D¹(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[4:6] atol=jatol rtol=jrtol
                @test D²(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[7:9] atol=jatol rtol=jrtol

                @test D³(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[10:12] atol=jatol rtol=jrtol

                # Velocity 
                @test D¹(t->ephem_vector6(ephj, cid, tid, t), tj) ≈ yj4[4:9] atol=jatol rtol=jrtol
                @test D²(t->ephem_vector6(ephj, cid, tid, t), tj) ≈ yj4[7:12]  atol=jatol rtol=jrtol

                # Acceleration 
                @test D¹(t->ephem_vector9(ephj, cid, tid, t), tj) ≈ yj4[4:12]  atol=jatol rtol=jrtol

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
    end

end;


