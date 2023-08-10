
using DiffRules
using ForwardDiff 

module TestDiffRules 

    using DiffRules 

    function test_diffdef(x, y, z)
        return x - 2y + 3z^2
    end 

    function d_test_diffdef(z)
        return cos(z)
    end
    
    using .TestDiffRules
    
    DiffRules.@define_diffrule TestDiffRules.test_diffdef(x, y, z) = :NaN, :NaN, :(d_test_diffdef($z))
    
end

# using .TestDiffRules

function test_fd(x, y, z)
    return x - 2y + 3z^2
end

""" THIS IS HOW YOU IMPLEMENT THE DERIVATION RULES FOR FORWARDDIFF"""
function test_fd(x, y, z::Dual{T}) where T 
    c = cos(value(z))
    Dual{T}(test_fd(x, y, value(z)), c*partials(z))
end

using BenchmarkTools

h = x -> test_fd(2, 3, x)
@benchmark ForwardDiff.derivative($h, π/2)
