
test_dir = artifact"testdata"
DJ2000 = 2451545

@testset "SPK Type 1" verbose=true begin 

    kernel = joinpath(test_dir, "spk1_ex1.bsp")

    # TODO: missing segment with epoch directory

    ephj = EphemerisProvider(kernel);
    furnsh(kernel)

    desc_id = rand(1:length(ephj.files[1].desc))

    desc = ephj.files[1].desc[desc_id]
    t1j, t2j = desc.tstart, desc.tend

    head = ephj.files[1].seglist.spk1[desc_id].head
    
    # Center and target bodies 
    cid = Int(desc.cid) 
    tid = Int(desc.tid)

    # Check errors 
    @test_throws jEphem.EphemerisError ephem_vector9(ephj, cid, tid, t1j)
    @test_throws jEphem.EphemerisError ephem_vector12(ephj, cid, tid, t1j)

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

        tj = j == 1 ? t1j : (j == 2 ? t2j : rand(ep))
        tc = tj/86400

        # Test with Julia
        yj1 = ephem_vector3(ephj, cid, tid, tj);
        yj2 = ephem_vector6(ephj, cid, tid, tj);

        # Test with SPICE
        ys1 = spkpos("$tid", tj, "J2000", "NONE", "$cid")[1]
        ys2 = spkezr("$tid", tj, "J2000", "NONE", "$cid")[1]

        @test yj1 ≈ ys1 atol=1e-13 rtol=1e-13
        @test yj2 ≈ ys2 atol=1e-13 rtol=1e-13

        # Test if AUTODIFF works 
        @test D¹(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj2[4:end] atol=1e-9 rtol=1e-9

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