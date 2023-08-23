
test_dir = artifact"testdata"
DJ2000 = 2451545

@testset "SPK Type 2" verbose=true begin 

    # --------------------------------------
    #           SPK files 
    # --------------------------------------

    kernel = joinpath(test_dir, "de421.bsp")

    ephj = EphemerisProvider(kernel);
    ephc = CalcephProvider(kernel);

    rec = ephem_spk_records(ephj)[2]
    t1j, t2j = rec.t_start[1], rec.t_end[1]

    # Test outside the time domain 
    @test_throws ephem_vector3(ephj, 399, 1,  t1j-1000) jEphem.EphemerisError 
    
    # Test values 
    yc1 = zeros(3);
    yc2 = zeros(6);
    yc3 = zeros(9);
    yc4 = zeros(12);

    ep = t1j:1:t2j
    for _ in 1:1000
        tj = rand(ep)
        tc = tj/86400

        yj1 = ephem_vector3(ephj, 399, 1,  tj);
        yj2 = ephem_vector6(ephj, 399, 1,  tj);
        yj3 = ephem_vector9(ephj, 399, 1,  tj);
        yj4 = ephem_vector12(ephj, 399, 1,  tj);
    
        jEphem.ephem_compute!(yc1, ephc, DJ2000, tc, 1,  399, 0);
        jEphem.ephem_compute!(yc2, ephc, DJ2000, tc, 1,  399, 1);
        jEphem.ephem_compute!(yc3, ephc, DJ2000, tc, 1,  399, 2);
        jEphem.ephem_compute!(yc4, ephc, DJ2000, tc, 1,  399, 3);
        
        # Test against CALCEPH
        @test yj1 ≈ yc1 atol=1e-9 rtol=1e-9
        @test yj2 ≈ yc2 atol=1e-9 rtol=1e-9
        @test yj3 ≈ yc3 atol=1e-9 rtol=1e-9
        @test yj4 ≈ yc4 atol=1e-9 rtol=1e-9

        # Test if AUTODIFF works 
        # Position
        @test D¹(t->ephem_vector3(ephj, 399, 1,  t), tj) ≈ yj4[4:6] atol=1e-9 rtol=1e-9
        @test D²(t->ephem_vector3(ephj, 399, 1,  t), tj) ≈ yj4[7:9] atol=1e-9 rtol=1e-9
        @test D³(t->ephem_vector3(ephj, 399, 1,  t), tj) ≈ yj4[10:12] atol=1e-9 rtol=1e-9

        # Velocity 
        @test D¹(t->ephem_vector6(ephj, 399, 1,  t), tj) ≈ yj4[4:9] atol=1e-9 rtol=1e-9
        @test D²(t->ephem_vector6(ephj, 399, 1,  t), tj) ≈ yj4[7:12] atol=1e-9 rtol=1e-9

        # Acceleration 
        @test D¹(t->ephem_vector9(ephj, 399, 1,  t), tj) ≈ yj4[4:12] atol=1e-9 rtol=1e-9

    end

    # --------------------------------------
    #           PCK files 
    # --------------------------------------

    kernel = joinpath(test_dir, "pa421.bsp")

    ephj = EphemerisProvider(kernel);
    ephc = CalcephProvider(kernel);

    rec = ephem_pck_records(ephj)[1]
    t1j, t2j = rec.t_start[1], rec.t_end[1]

    # Test outside the time domain 
    @test_throws ephem_rotation3(ephj, 31006, 1, t1j-1000) jEphem.EphemerisError 
    
    # Test values 
    yc1 = zeros(3);
    yc2 = zeros(6);
    yc3 = zeros(9);
    yc4 = zeros(12);

    ep = t1j:1:t2j
    for _ in 1:1000
        tj = rand(ep)
        tc = tj/86400

        yj1 = ephem_rotation3(ephj, 31006, 3, tj);
        yj2 = ephem_rotation6(ephj, 31006, 3, tj);
        yj3 = ephem_rotation9(ephj, 31006, 3, tj);
        yj4 = ephem_rotation12(ephj, 31006, 3, tj);
    
        jEphem.ephem_orient!(yc1, ephc, DJ2000, tc, 31006, 1, 0);
        jEphem.ephem_orient!(yc2, ephc, DJ2000, tc, 31006, 1, 1);
        jEphem.ephem_orient!(yc3, ephc, DJ2000, tc, 31006, 1, 2);
        jEphem.ephem_orient!(yc4, ephc, DJ2000, tc, 31006, 1, 3);
        
        # Test against CALCEPH
        @test yj1 ≈ yc1 atol=1e-9 rtol=1e-9
        @test yj2 ≈ yc2 atol=1e-9 rtol=1e-9
        @test yj3 ≈ yc3 atol=1e-9 rtol=1e-9
        @test yj4 ≈ yc4 atol=1e-9 rtol=1e-9

        # Test if AUTODIFF works 
        # Position
        @test D¹(t->ephem_rotation3(ephj, 31006, 1, t), tj) ≈ yj4[4:6] atol=1e-9 rtol=1e-9
        @test D²(t->ephem_rotation3(ephj, 31006, 1, t), tj) ≈ yj4[7:9] atol=1e-9 rtol=1e-9
        @test D³(t->ephem_rotation3(ephj, 31006, 1, t), tj) ≈ yj4[10:12] atol=1e-9 rtol=1e-9

        # Velocity 
        @test D¹(t->ephem_rotation6(ephj, 31006, 1, t), tj) ≈ yj4[4:9] atol=1e-9 rtol=1e-9
        @test D²(t->ephem_rotation6(ephj, 31006, 1, t), tj) ≈ yj4[7:12] atol=1e-9 rtol=1e-9

        # Acceleration 
        @test D¹(t->ephem_rotation9(ephj, 31006, 1, t), tj) ≈ yj4[4:12] atol=1e-9 rtol=1e-9

    end
end



