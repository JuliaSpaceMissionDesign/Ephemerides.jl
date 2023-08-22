
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

    y = zeros(6)
    yc = zeros(6);

    # Point ephemeris tests
    @test_throws jEphem.EphemerisError jEphem.ephem_compute!(y, de421, 0.0, 0.0, 301, 399, 1)
    jEphem.ephem_compute!(y, de421, DJ2000, 0.0, 301, 3, 1)
    
    jEphem.ephem_compute!(yc, ephc, DJ2000, 0.0, 301, 3, 1)
    @test yc ≈ y atol=1e-11 rtol=1e-11

    # # Axes ephemeris tests 
    # @test_throws jEphem.EphemerisError jEphem.ephem_orient!(y, pa421, 0.0, 0.0, 301, 1)
    # jEphem.ephem_orient!(y, pa421, DJ2000, 0.0, 31006, 1)

    # jEphem.ephem_orient!(yc, epho, DJ2000, 0.0, 31006, 1)
    # @test yc ≈ y atol=1e-11 rtol=1e-11

end;