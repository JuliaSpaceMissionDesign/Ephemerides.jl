
@testset "TwoBodyUtils" verbose=true begin 

    deg = 15 
    
    p = zeros(2(deg-1))
    @inbounds for j in eachindex(p)
        p[j] = 1/(j*(j+1))
    end

    function stumpff_series(x, p)

        # Manually compute the Maclaurin series 
        n = length(p)

        # Compute the series expansion for C₃:
        c3 = one(x) 
        @inbounds for j = n:-2:4 
            c3 = 1 - x*p[j]*c3
        end
        c3 *= p[2]

        # Compute the series expansion for C₂: 
        c2 = one(x)
        @inbounds for j = n-1:-2:3 
            c2 = 1 - x*p[j]*c2 
        end 
        c2 *= p[1]

        # Compute C₀ and C₁ from the recursion formula: 
        c1 = 1 - x*c3 
        c0 = 1 - x*c2

        return c0, c1, c2, c3

    end

    x = -10:1e-3:10

    for j in eachindex(x)

        a1, a2, a3, a4 = stumpff_series(BigFloat(x[j]), p)
        b1, b2, b3, b4 = Ephemerides.stumpff(x[j], p)

        @test a1 ≈ b1 atol=1e-15 rtol=1e-15
        @test a2 ≈ b2 atol=1e-15 rtol=1e-15
        @test a3 ≈ b3 atol=1e-15 rtol=1e-15
        @test a4 ≈ b4 atol=1e-15 rtol=1e-15

    end

end;