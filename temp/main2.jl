using ForwardDiff
using ForwardDiff: derivative
using PreallocationTools

using BenchmarkTools

randmat = rand(5, 3)
sto = similar(randmat)
stod = DiffCache(sto)

@views function claytonsample!(cache, τ, α, randmat)
    
    sto = get_tmp(cache, τ)
    sto .= randmat
    τ == 0 && return sto

    print(typeof(sto))

    n = size(sto, 1)
    for i in 1:n
        v = sto[i, 2]
        u = sto[i, 1]
        sto[i, 1] = (1 - u^(-τ) + u^(-τ)*v^(-(τ/(1 + τ))))^(-1/τ)*α
        sto[i, 2] = (1 - u^(-τ) + u^(-τ)*v^(-(τ/(1 + τ))))^(-1/τ)
    end

    return sto

end

ForwardDiff.derivative(τ -> claytonsample!(stod, τ, 0.0, randmat), 0.3)

@benchmark ForwardDiff.derivative(τ -> claytonsample!($stod, τ, $(0.0), $randmat), 0.3)

struct TestCache 
    x::Vector{Float64}
end

cac = TestCache(zeros(10));

function get_cache(cache::TestCache)
    return @views cache.x
end

function fcn(cac)

    x = get_cache(cac)
    x[1] = 100
    nothing 

end


@benchmark fcn($cac)

u = zeros(3)

test = DiffCache(u)
length(test.dual_du)

ForwardDiff.pickchunksize(length(u))


struct TestDualCache 
    u::Vector{Float64}
end

@inline get_u(cache::TestDualCache, t) = @views cache.u
function get_u(cache::TestDualCache, t::T) where {T <: ForwardDiff.Dual}
    reinterpret(T, @views(cache.u[:]))
end

u = Float64[1, 2, 3, 4, 5];
tcache = TestDualCache([1, 2, 3, 4, 5])
dcache = DiffCache(zeros(5))

function fcn(cache::TestDualCache, t)

    u = get_u(cache, t)
    u[1] = 4

    return u
end



function fcntest(cache::DiffCache, t::T) where {T}

    u = get_tmp(cache, t)
    u[1] = 1.0 

    for i = 2:lastindex(u)
        u[i] = t*u[i-1]
    end

    x, y, z = T(0.0), T(0.0), T(0.0)
    for i in eachindex(u)
        x += u[i]
        y += 2u[i]
        z += 3u[i]
    end

    return SA[x, y, z]
end

@inline D¹(f, t) = derivative(f, t) 
@inline D²(f, t) = derivative(τ->derivative(f, τ), t)
@inline D³(f, t) = derivative(κ->derivative(τ->derivative(f, τ), κ), t)

fcntest(dcache, 2.0)
D¹(t->fcntest(dcache, t), 2.0)
D²(t->fcntest(dcache, t), 2.0)
D³(t->fcntest(dcache, t), 2.0)

@benchmark fcntest($dcache, $(2.0))
@benchmark D¹(t->fcntest($dcache, t), $(2.0))
@benchmark D²(t->fcntest($dcache, t), $(2.0))
@benchmark D³(t->fcntest($dcache, t), $(2.0))
