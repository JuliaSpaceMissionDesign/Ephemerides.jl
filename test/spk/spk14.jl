
test_dir = artifact"testdata"
DJ2000 = 2451545

@testset "SPK Type 14" verbose=true begin 

    # The first kernel has a directory epoch, the second does not
    kernel = joinpath(test_dir, "spk14_ex1.bsp")

    ephj = EphemerisProvider(kernel);
    furnsh(kernel)

    cid = 10 
    tids = [2000, 2001, 2002]

    for tid in tids 
        
        link = ephj.spklinks[cid][tid][1]
        
        head = ephj.files[link.fid].seglist.spk14[link.eid].head 

        desc = ephj.spklinks[cid][tid][1].desc 
        t1j, t2j = desc.tstart, desc.tend

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
                tj = min(t2j, max(rand(head.epochs), t1j))
            elseif j < 500
                # Test directory handling close to the borders
                tj = min(t2j, max(rand(head.epochs) + randn(), t1j))
            else             
                tj = t1j + (t2j-t1j)*rand()
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
            # Position (position doesn't work exactly because they have different coefficients)
            # @test D¹(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[4:6] atol=1e-7 rtol=1e-12
            # @test D²(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[7:9] atol=1e-7 rtol=1e-12
            # @test D³(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[10:12] atol=1e-7 rtol=1e-12

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

    end

    kclear()

end;
