
test_dir = artifact"testdata"
DJ2000 = 2451545.0

@testset "JSMD Interfaces" verbose=true begin 

    path_de421 = joinpath(test_dir, "de421.bsp")
    path_pa421 = joinpath(test_dir, "pa421.bpc")

    de421 = EphemerisProvider(path_de421)
    pa421 = jEphem.load(EphemerisProvider, path_pa421)

    kern = jEphem.load(EphemerisProvider, [path_pa421, path_de421])

    points = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 199, 299, 301, 399, 499]
    axes = [1, 31006]

    # CALCEPH.jl files 
    ephc = CalcephProvider(path_de421)
    epho = CalcephProvider(path_pa421)

    # Time properties tests
    @test jEphem.ephem_timescale(de421) == 1
    @test jEphem.ephem_timescale(pa421) == 1

    @test jEphem.ephem_timespan(de421) == jEphem.ephem_timespan(ephc)
    @test jEphem.ephem_timespan(pa421) == jEphem.ephem_timespan(epho)

    @test jEphem.ephem_available_axes(de421) == Int64[]
    @test jEphem.ephem_available_points(pa421) == Int64[]

    @test sort(jEphem.ephem_available_axes(pa421)) == axes
    @test sort(jEphem.ephem_available_points(de421)) == points

    @test sort(jEphem.ephem_available_axes(kern)) == axes
    @test sort(jEphem.ephem_available_points(kern)) == points

    # Data records tests
    @test jEphem.ephem_position_records(pa421) == jEphem.EphemPointRecord[]
    @test jEphem.ephem_orient_records(de421) == jEphem.EphemAxesRecord[]

    prec = jEphem.ephem_position_records(de421)
    crec = sort(jEphem.ephem_position_records(ephc), by = x -> x.target, rev = true)

    for (rec, recc) in zip(prec, crec) 
        
        @test ((rec.center == recc.center && rec.target == recc.target) || 
               (rec.center == recc.target && rec.target == recc.center)) 

        @test rec.axes == recc.axes 
        @test rec.jd_start == recc.jd_start 
        @test rec.jd_stop == recc.jd_stop
    end

    orec = jEphem.ephem_orient_records(pa421)
    corec = sort(jEphem.ephem_orient_records(epho), rev = true)

    for (rec, recc) in zip(orec, corec) 
        @test ((rec.target == recc.target && rec.axes == recc.axes) || 
               (rec.axes == recc.target && rec.target == recc.axes))

        @test rec.jd_start == recc.jd_start 
        @test rec.jd_stop == recc.jd_stop
    end

    yj1, yc1 = zeros(3), zeros(3)
    yj2, yc2 = zeros(6), zeros(6)
    yj3, yc3 = zeros(9), zeros(9)
    yj4, yc4 = zeros(12), zeros(12)

    # Point ephemeris tests
    @test_throws jEphem.EphemerisError jEphem.ephem_compute!(yj1, de421, 0.0, 0.0, 301, 399, 0)

    jEphem.ephem_compute!(yj1, de421, DJ2000, 0.0, 301, 3, 0)
    jEphem.ephem_compute!(yc1, ephc, DJ2000, 0.0, 301, 3, 0)

    jEphem.ephem_compute!(yj2, de421, DJ2000, 0.0, 301, 3, 1)
    jEphem.ephem_compute!(yc2, ephc, DJ2000, 0.0, 301, 3, 1)

    jEphem.ephem_compute!(yj3, de421, DJ2000, 0.0, 301, 3, 2)
    jEphem.ephem_compute!(yc3, ephc, DJ2000, 0.0, 301, 3, 2)

    jEphem.ephem_compute!(yj4, de421, DJ2000, 0.0, 301, 3, 3)
    jEphem.ephem_compute!(yc4, ephc, DJ2000, 0.0, 301, 3, 3)

    @test yc1 ≈ yj1 atol=1e-11 rtol=1e-11
    @test yc2 ≈ yj2 atol=1e-11 rtol=1e-11
    @test yc3 ≈ yj3 atol=1e-11 rtol=1e-11
    @test yc4 ≈ yj4 atol=1e-11 rtol=1e-11

    # Axes ephemeris tests 
    @test_throws jEphem.EphemerisError jEphem.ephem_orient!(yj1, pa421, 0.0, 0.0, 301, 1, 0)

    jEphem.ephem_orient!(yj1, pa421, DJ2000, 0.0, 31006, 1, 0)
    jEphem.ephem_orient!(yc1, epho, DJ2000, 0.0, 31006, 1, 0)

    jEphem.ephem_orient!(yj2, pa421, DJ2000, 0.0, 31006, 1, 1)
    jEphem.ephem_orient!(yc2, epho, DJ2000, 0.0, 31006, 1, 1)

    jEphem.ephem_orient!(yj3, pa421, DJ2000, 0.0, 31006, 1, 2)
    jEphem.ephem_orient!(yc3, epho, DJ2000, 0.0, 31006, 1, 2)

    jEphem.ephem_orient!(yj4, pa421, DJ2000, 0.0, 31006, 1, 3)
    jEphem.ephem_orient!(yc4, epho, DJ2000, 0.0, 31006, 1, 3)

    @test yc1 ≈ yj1 atol=1e-11 rtol=1e-11
    @test yc2 ≈ yj2 atol=1e-11 rtol=1e-11
    @test yc3 ≈ yj3 atol=1e-11 rtol=1e-11
    @test yc4 ≈ yj4 atol=1e-11 rtol=1e-11

end;