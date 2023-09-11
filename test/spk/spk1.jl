
test_dir = artifact"testdata"
DJ2000 = 2451545

@testset "SPK Type 1" verbose=true begin 

    kernel = joinpath(test_dir, "example1spk_seg1.bsp")

    ephj = EphemerisProvider(kernel);
    furnsh(kernel)

    desc = rand(ephj.files[1].desc)
    t1j, t2j = desc.tstart, desc.tend
    
    # Center and target bodies 
    cid = Int(desc.cid) 
    tid = Int(desc.tid)

    ep = t1j:1:t2j
    for j in 1:2000

        tj = j == 1 ? t1j : (j == 2 ? t2j : rand(ep))
        tc = tj/86400

        # Test with Julia
        yj1 = ephem_vector3(ephj, cid, tid, tj);
        yj2 = ephem_vector6(ephj, cid, tid, tj);

        # Test with SPICE
        ys1 = spkpos("$tid", tj, "J2000", "NONE", "$cid")[1]
        ys2 = spkezr("$tid", tj, "J2000", "NONE", "$cid")[1]

        @test yj1 ≈ ys1 atol=1e-9 rtol=1e-13
        @test yj2 ≈ ys2 atol=1e-9 rtol=1e-13

        # Test if AUTODIFF works 
        @test D¹(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj2[4:end] atol=1e-9 rtol=1e-9

        # TODO: implement the acceleration and jerk. This functions below cannot be tested!
        # D²(t->ephem_vector3(ephj, cid, tid, t), tj)
        # D³(t->ephem_vector3(ephj, cid, tid, t), tj)

        # D¹(t->ephem_vector6(ephj, cid, tid, t), tj)
        # D²(t->ephem_vector6(ephj, cid, tid, t), tj)
        # D³(t->ephem_vector6(ephj, cid, tid, t), tj)

    end

    kclear()

end