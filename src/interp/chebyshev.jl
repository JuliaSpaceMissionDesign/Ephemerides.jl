
"""
    chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int)

Evaluate a sum of Cheybyshev polynomials of the first kind at `t` using a 
recursive algorithm. It simultenously evalutes the 3 state components. `idx` is the 
index of the starting row (in 0-based notation) in the matrix of coefficients `cₖ`.

!!! note 
    `x` is a re-work of the actual ascissa value that lies between [-1, 1]

### See Also 
See also [`∂chebyshev`](@ref), [`∂²chebyshev`](@ref) and [`∂³chebyshev`](@ref)

"""
function chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int)

    # Retrieve the work buffer 
    Tₙ = get_buffer(cache, 1, t)

    # idx is 0 or 3
    ix = 1 + idx 
    iy = 2 + idx 
    iz = 3 + idx 

    @inbounds begin 

        Tₙ[1] = 1
        Tₙ[2] = t 
        Tₙ[3] = 2t*t - 1

        x = cₖ[ix, 1] + t*cₖ[ix, 2] + Tₙ[3]*cₖ[ix, 3]
        y = cₖ[iy, 1] + t*cₖ[iy, 2] + Tₙ[3]*cₖ[iy, 3]
        z = cₖ[iz, 1] + t*cₖ[iz, 2] + Tₙ[3]*cₖ[iz, 3]
        
        for j = 4:N
            
            Tₙ[j] = 2*t*Tₙ[j-1] - Tₙ[j-2]
        
            x += Tₙ[j]*cₖ[ix, j]
            y += Tₙ[j]*cₖ[iy, j]
            z += Tₙ[j]*cₖ[iz, j]
        
        end
        
    end

    return x, y, z

end 


"""
    ∂chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int, Δt)

Evaluate a sum of Cheybyshev polynomials of the first kind and its derivative at `t` 
using a recursive algorithm. It simultenously evalutes the 3 state components. `idx` 
is the index of the starting row (in 0-based notation) in the matrix of coefficients `cₖ`.

!!! note 
    `x` is a re-work of the actual ascissa value that lies between [-1, 1]

### See Also 
See also [`chebyshev`](@ref), [`∂²chebyshev`](@ref) and [`∂³chebyshev`](@ref)

"""
function ∂chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int, Δt::Number)

    # Retrieve the work buffer 
    Tₙ  = get_buffer(cache, 1, t)
    dTₙ = get_buffer(cache, 2, t)

    # idx is 0 or 3
    ix = 1 + idx 
    iy = 2 + idx 
    iz = 3 + idx 

    @inbounds begin 

        Tₙ[1] = 1
        Tₙ[2] = t 
        Tₙ[3] = 2t*t - 1

        dTₙ[2] = 0 
        dTₙ[2] = 1 
        dTₙ[3] = 4t

        x = cₖ[ix, 1] + t*cₖ[ix, 2] + Tₙ[3]*cₖ[ix, 3]
        y = cₖ[iy, 1] + t*cₖ[iy, 2] + Tₙ[3]*cₖ[iy, 3]
        z = cₖ[iz, 1] + t*cₖ[iz, 2] + Tₙ[3]*cₖ[iz, 3]

        vx = cₖ[ix, 2] + dTₙ[3]*cₖ[ix, 3]
        vy = cₖ[iy, 2] + dTₙ[3]*cₖ[iy, 3]
        vz = cₖ[iz, 2] + dTₙ[3]*cₖ[iz, 3]
        
        for j = 4:N

            Tₙ[j]  = 2*t*Tₙ[j-1] - Tₙ[j-2]
            dTₙ[j] = 2*t*dTₙ[j-1] + 2*Tₙ[j-1] - dTₙ[j-2]
        
            x += Tₙ[j]*cₖ[ix, j]
            y += Tₙ[j]*cₖ[iy, j]
            z += Tₙ[j]*cₖ[iz, j]

            vx += dTₙ[j]*cₖ[ix, j]
            vy += dTₙ[j]*cₖ[iy, j]
            vz += dTₙ[j]*cₖ[iz, j]
        end
        
    end
    
    vx *= Δt 
    vy *= Δt 
    vz *= Δt

    return x, y, z, vx, vy, vz

end 


"""
    ∂²chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int, Δt)

Evaluate a sum of Cheybyshev polynomials of the first kind and its two derivatives at `t` 
using a recursive algorithm. It simultenously evalutes the 3 state components. `idx` 
is the index of the starting row (in 0-based notation) in the matrix of coefficients `cₖ`.

!!! note 
    `x` is a re-work of the actual ascissa value that lies between [-1, 1]

### See Also 
See also [`chebyshev`](@ref), [`∂chebyshev`](@ref) and [`∂³chebyshev`](@ref)

"""
function ∂²chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int, Δt)

    # Retrieve the work buffer 
    Tₙ   = get_buffer(cache, 1, t)
    dTₙ  = get_buffer(cache, 2, t)
    ddTₙ = get_buffer(cache, 3, t)

    # idx is 0 or 3
    ix = 1 + idx 
    iy = 2 + idx 
    iz = 3 + idx 

    @inbounds begin 

        Tₙ[1] = 1
        Tₙ[2] = t 
        Tₙ[3] = 2t*t - 1

        dTₙ[1] = 0 
        dTₙ[2] = 1 
        dTₙ[3] = 4t

        ddTₙ[1] = 0
        ddTₙ[2] = 0
        ddTₙ[3] = 4

        x = cₖ[ix, 1] + t*cₖ[ix, 2] + Tₙ[3]*cₖ[ix, 3]
        y = cₖ[iy, 1] + t*cₖ[iy, 2] + Tₙ[3]*cₖ[iy, 3]
        z = cₖ[iz, 1] + t*cₖ[iz, 2] + Tₙ[3]*cₖ[iz, 3]

        vx = cₖ[ix, 2] + dTₙ[3]*cₖ[ix, 3]
        vy = cₖ[iy, 2] + dTₙ[3]*cₖ[iy, 3]
        vz = cₖ[iz, 2] + dTₙ[3]*cₖ[iz, 3]

        ax = 4cₖ[ix, 3]
        ay = 4cₖ[iy, 3]
        az = 4cₖ[iz, 3]
        
        for j = 4:N

            Tₙ[j]   = 2*t*Tₙ[j-1]                - Tₙ[j-2]
            dTₙ[j]  = 2*t*dTₙ[j-1]  + 2*Tₙ[j-1]  - dTₙ[j-2]
            ddTₙ[j] = 2*t*ddTₙ[j-1] + 4*dTₙ[j-1] - ddTₙ[j-2]
            
            x += Tₙ[j]*cₖ[ix, j]
            y += Tₙ[j]*cₖ[iy, j]
            z += Tₙ[j]*cₖ[iz, j]

            vx += dTₙ[j]*cₖ[ix, j]
            vy += dTₙ[j]*cₖ[iy, j]
            vz += dTₙ[j]*cₖ[iz, j]

            ax += ddTₙ[j]*cₖ[ix, j]
            ay += ddTₙ[j]*cₖ[iy, j]
            az += ddTₙ[j]*cₖ[iz, j]

        end
        
    end

    vx *= Δt
    vy *= Δt
    vz *= Δt

    Δt² = Δt*Δt

    ax *= Δt²
    ay *= Δt²
    az *= Δt²

    return x, y, z, vx, vy, vz, ax, ay, az

end 

"""
    ∂³chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int, Δt)

Evaluate a sum of Cheybyshev polynomials of the first kind and its three derivatives at `t` 
using a recursive algorithm. It simultenously evalutes the 3 state components. `idx` 
is the index of the starting row (in 0-based notation) in the matrix of coefficients `cₖ`.

!!! note 
    `x` is a re-work of the actual ascissa value that lies between [-1, 1]

### See Also 
See also [`chebyshev`](@ref), [`∂chebyshev`](@ref) and [`∂²chebyshev`](@ref)

"""
function ∂³chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int, Δt)

    # Retrieve the work buffer 
    Tₙ    = get_buffer(cache, 1, t)
    dTₙ   = get_buffer(cache, 2, t)
    ddTₙ  = get_buffer(cache, 3, t)
    dddTₙ = get_buffer(cache, 4, t)

    # idx is 0 or 3
    ix = 1 + idx 
    iy = 2 + idx 
    iz = 3 + idx 

    @inbounds begin 

        Tₙ[1] = 1
        Tₙ[2] = t 
        Tₙ[3] = 2t*t - 1

        dTₙ[1] = 0 
        dTₙ[2] = 1 
        dTₙ[3] = 4t

        ddTₙ[1] = 0
        ddTₙ[2] = 0
        ddTₙ[3] = 4

        dddTₙ[1] = 0
        dddTₙ[2] = 0
        dddTₙ[3] = 0

        x = cₖ[ix, 1] + t*cₖ[ix, 2] + Tₙ[3]*cₖ[ix, 3]
        y = cₖ[iy, 1] + t*cₖ[iy, 2] + Tₙ[3]*cₖ[iy, 3]
        z = cₖ[iz, 1] + t*cₖ[iz, 2] + Tₙ[3]*cₖ[iz, 3]

        vx = cₖ[ix, 2] + dTₙ[3]*cₖ[ix, 3]
        vy = cₖ[iy, 2] + dTₙ[3]*cₖ[iy, 3]
        vz = cₖ[iz, 2] + dTₙ[3]*cₖ[iz, 3]

        ax = 4cₖ[ix, 3]
        ay = 4cₖ[iy, 3]
        az = 4cₖ[iz, 3]
        
        jx = dddTₙ[1]
        jy = dddTₙ[1]
        jz = dddTₙ[1]

        for j = 4:N

            Tₙ[j]    = 2t*Tₙ[j-1]                  - Tₙ[j-2]
            dTₙ[j]   = 2t*dTₙ[j-1]   + 2*Tₙ[j-1]   - dTₙ[j-2]
            ddTₙ[j]  = 2t*ddTₙ[j-1]  + 4*dTₙ[j-1]  - ddTₙ[j-2]
            dddTₙ[j] = 2t*dddTₙ[j-1] + 6*ddTₙ[j-1] - dddTₙ[j-2]
            
            x += Tₙ[j]*cₖ[ix, j]
            y += Tₙ[j]*cₖ[iy, j]
            z += Tₙ[j]*cₖ[iz, j]

            vx += dTₙ[j]*cₖ[ix, j]
            vy += dTₙ[j]*cₖ[iy, j]
            vz += dTₙ[j]*cₖ[iz, j]

            ax += ddTₙ[j]*cₖ[ix, j]
            ay += ddTₙ[j]*cₖ[iy, j]
            az += ddTₙ[j]*cₖ[iz, j]

            jx += dddTₙ[j]*cₖ[ix, j]
            jy += dddTₙ[j]*cₖ[iy, j]
            jz += dddTₙ[j]*cₖ[iz, j]
        
        end
        
    end

    vx *= Δt
    vy *= Δt
    vz *= Δt

    Δt² = Δt*Δt

    ax *= Δt²
    ay *= Δt²
    az *= Δt²

    Δt³ = Δt²*Δt

    jx *= Δt³
    jy *= Δt³
    jz *= Δt³

    return x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz

end 