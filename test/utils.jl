
function test_vsep(x, y)

    nx, ny = norm(x), norm(y)

    ux = x/nx
    uy = y/ny

    if dot(ux, uy) > 0 
        l = norm(ux - uy)
        return 2asin(l/2)
    elseif dot(ux, uy) < 0 
        l = norm(ux + uy)
        return π - 2asin(l/2)
    else 
        return π/2
    end

end

@testset "Vector Utils" verbose=true begin 

    for _ in 1:1000
        x = rand(3) 
        y = rand(3)

        @test norm(x) ≈ Ephemerides.vnorm(x)           atol=1e-13 rtol=1e-13
        @test x/norm(x) ≈ Ephemerides.vhat(x)          atol=1e-13 rtol=1e-13
        @test dot(x, y) ≈ Ephemerides.vdot(x, y)       atol=1e-13 rtol=1e-13
        @test cross(x, y) ≈ Ephemerides.vcross(x, y)   atol=1e-13 rtol=1e-13
        @test vsep(x, y) ≈ Ephemerides.vsep(x, y)      atol=1e-13 rtol=1e-13
    end

    k = [0, 0, 1]
    @test Ephemerides.vrot([1, 2, 3], k, π/2) ≈ [-2, 1, 3] atol=1e-13 rtol=1e-13
    @test Ephemerides.vrot([1, 0, 0], k, π/2) ≈ [0, 1, 0] atol=1e-13 rtol=1e-13
    @test Ephemerides.vrot([0, 1, 0], k, π/2) ≈ [-1, 0, 0] atol=1e-13 rtol=1e-13

    k = [0, 1, 0]
    @test Ephemerides.vrot([1, 2, 3], k, π/2) ≈ [3, 2, -1] atol=1e-13 rtol=1e-13
    @test Ephemerides.vrot([1, 0, 0], k, π/2) ≈ [0, 0, -1] atol=1e-13 rtol=1e-13
    @test Ephemerides.vrot([0, 1, 0], k, π/2) ≈ [0, 1, 0] atol=1e-13 rtol=1e-13

end;