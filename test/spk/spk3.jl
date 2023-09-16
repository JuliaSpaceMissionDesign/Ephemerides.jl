
test_dir = artifact"testdata"
DJ2000 = 2451545

@testset "SPK Type 3" verbose=true begin 

    kernel = joinpath(test_dir, "example1spk_seg3.bsp")

    ephj = EphemerisProvider(kernel);
    ephc = CalcephProvider(kernel);
    furnsh(kernel)

    desc = rand(ephj.files[1].desc)
    t1j, t2j = desc.tstart, desc.tend

    cid = Int(desc.cid)
    tid = Int(desc.tid)

    # Test values 
    yc1 = zeros(3);
    yc2 = zeros(6);
    yc3 = zeros(9);
    yc4 = zeros(12);

    ep = t1j:1:t2j
    for j in 1:2000

        tj = j == 1 ? t1j : (j == 2 ? t2j : rand(ep))
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
        @test yj1 ≈ yc1 atol=1e-9 rtol=1e-11
        @test yj2 ≈ yc2 atol=1e-9 rtol=1e-11
        @test yj3 ≈ yc3 atol=1e-9 rtol=1e-11
        @test yj4 ≈ yc4 atol=1e-9 rtol=1e-11

        # Test against SPICE 
        @test ys1 ≈ yj1 atol=1e-13 rtol=1e-14
        @test ys2 ≈ yj2 atol=1e-13 rtol=1e-14
        
        # Test if AUTODIFF works 
        # Position (position doesn't work exactly because they have different coefficients)
        # @test D¹(t->ephem_vector3(ephj, cid, tid,  t), tj) ≈ yj4[4:6] atol=1e-5 rtol=1e-5
        # @test D²(t->ephem_vector3(ephj, cid, tid,  t), tj) ≈ yj4[7:9] atol=1e-5 rtol=1e-5
        # @test D³(t->ephem_vector3(ephj, cid, tid,  t), tj) ≈ yj4[10:12] atol=1e-5 rtol=1e-5

        # Velocity 
        @test D¹(t->ephem_vector6(ephj, cid, tid,  t), tj)[4:6] ≈ yj4[7:9] atol=1e-9 rtol=1e-12
        @test D²(t->ephem_vector6(ephj, cid, tid,  t), tj)[4:6] ≈ yj4[10:12] atol=1e-8 rtol=1e-12

        # Acceleration 
        @test D¹(t->ephem_vector9(ephj, cid, tid,  t), tj)[4:9] ≈ yj4[7:12] atol=1e-8 rtol=1e-12

    end

    kclear()

end



