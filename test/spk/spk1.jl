
test_dir = artifact"testdata"
DJ2000 = 2451545

@testset "SPK Type 1" verbose=true begin 

    kernel = joinpath(test_dir, "example1spk_seg1.bsp")

    ephj = EphemerisProvider(kernel);
    ephc = CalcephProvider(kernel);

    # Center and target bodies 
    cid = 0 
    tid = 2000001

    t1j, t2j, tcj = ephem_spk_timespan(ephj)

    # Check that the timespan is correct 
    t1c, t2c, tcc = jEphem.ephem_timespan(ephc)

    @test t1c == t1j/86400 + DJ2000 
    @test t2c == t2j/86400 + DJ2000
    @test tcc == tcj

    # Test values 
    yc1 = zeros(3)
    yc2 = zeros(6)

    ep = tj1:1:t2j
    for _ in 1:1000

        tj = rand(ep)
        tc = tj/86400

        # Test with Julia
        yj1 = ephem_vector3(ephj, cid, tid, tj);
        yj2 = ephem_vector6(ephj, cid, tid, tj);

        # Test with Calceph
        jEphem.ephem_compute!(yc1, ephc, DJ2000, tc, tid, cid, 0);
        jEphem.ephem_compute!(yc2, ephc, DJ2000, tc, tid, cid, 1);

        @test yj1 ≈ yc1 atol=1e-9 rtol=1e-9
        @test yj2 ≈ yc2 atol=1e-9 rtol=1e-9

        # Test if AUTODIFF works 
        @test D¹(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj2[4:end] atol=1e-9 rtol=1e-9

        # TODO: implement the acceleration and jerk. This functions below cannot be tested!
        # D²(t->ephem_vector3(ephj, cid, tid, t), tj)
        # D³(t->ephem_vector3(ephj, cid, tid, t), tj)

        # D¹(t->ephem_vector6(ephj, cid, tid, t), tj)
        # D²(t->ephem_vector6(ephj, cid, tid, t), tj)
        # D³(t->ephem_vector6(ephj, cid, tid, t), tj)

    end

end