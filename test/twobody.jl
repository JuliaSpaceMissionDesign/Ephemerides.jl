
@testset "TwoBodyUtils" verbose=true begin 

    deg = 15 
    
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

    function stumpff_series()
    end

end;