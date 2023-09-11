
test_dir = artifact"testdata"

DJ2000 = 2451545

@testset "SPK Type 18" verbose=true begin 

    kernel = joinpath(test_dir, "example1spk_seg18.bsp") # subtype = 1
    
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
                
        ys1 = spkpos("$tid", tj, "J2000", "NONE", "$cid")[1]
        ys2 = spkezr("$tid", tj, "J2000", "NONE", "$cid")[1]

        # Test against SPICE
        @test yj1 ≈ ys1 atol=1e-9 rtol=1e-14
        @test yj2 ≈ ys2 atol=1e-9 rtol=1e-14

        # Test if AUTODIFF works 
        # Position (doesn't work exactly because pos and vel have different set of coefficients)
        # @test D¹(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[4:6] atol=1e-9 rtol=1e-13
        # @test D²(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[7:9] atol=1e-9 rtol=1e-13
        # @test D³(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj4[10:12] atol=1e-9 rtol=1e-13

        # Velocity 
        @test D¹(t->ephem_vector6(ephj, cid, tid, t), tj)[4:6] ≈ yj3[7:9] atol=1e-9 rtol=1e-14
        @test D²(t->ephem_vector6(ephj, cid, tid, t), tj)[4:6] ≈ yj4[10:12] atol=1e-9 rtol=1e-14

        # Acceleration 
        @test D¹(t->ephem_vector9(ephj, cid, tid, t), tj)[4:9] ≈ yj4[7:12] atol=1e-9 rtol=1e-14

    end

    kclear()

end;


