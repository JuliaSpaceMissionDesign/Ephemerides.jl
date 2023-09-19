
"""
    chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int, ibuff::Int)

Evaluate a sum of Cheybyshev polynomials of the first kind at `t` using a 
recursive algorithm. It simultenously evalutes the 3 state components. `idx` is the 
index of the starting row (in 0-based notation) in the matrix of coefficients `cₖ`.
`ibuff` is the index of the first free buffer. 

!!! note 
    `x` is a re-work of the actual ascissa value that lies between [-1, 1]

### See Also 
See also [`∂chebyshev`](@ref), [`∂²chebyshev`](@ref) and [`∂³chebyshev`](@ref)

"""
function chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int, ibuff=1)

    # Retrieve the work buffer 
    Tₙ = get_buffer(cache, ibuff, t)

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
    ∂chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int, Δt, ibuff=1)

Evaluate a sum of Cheybyshev polynomials of the first kind and its derivative at `t` 
using a recursive algorithm. It simultenously evalutes the 3 state components. `idx` 
is the index of the starting row (in 0-based notation) in the matrix of coefficients `cₖ`.
`ibuff` is the index of the first free buffer. 

!!! note 
    `x` is a re-work of the actual ascissa value that lies between [-1, 1]

### See Also 
See also [`chebyshev`](@ref), [`∂²chebyshev`](@ref) and [`∂³chebyshev`](@ref)

"""
function ∂chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int, Δt::Number, ibuff=1)

    # Retrieve the work buffer 
    Tₙ  = get_buffer(cache, ibuff, t)
    dTₙ = get_buffer(cache, ibuff + 1, t)

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

    return x, y, z, Δt*vx, Δt*vy, Δt*vz

end 


"""
    ∂²chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int, Δt, ibuff=1)

Evaluate a sum of Cheybyshev polynomials of the first kind and its two derivatives at `t` 
using a recursive algorithm. It simultenously evalutes the 3 state components. `idx` 
is the index of the starting row (in 0-based notation) in the matrix of coefficients `cₖ`.
`ibuff` is the index of the first free buffer. 

!!! note 
    `x` is a re-work of the actual ascissa value that lies between [-1, 1]

### See Also 
See also [`chebyshev`](@ref), [`∂chebyshev`](@ref) and [`∂³chebyshev`](@ref)

"""
function ∂²chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int, Δt, ibuff=1)

    # Retrieve the work buffer 
    Tₙ   = get_buffer(cache, ibuff, t)
    dTₙ  = get_buffer(cache, ibuff + 1, t)
    ddTₙ = get_buffer(cache, ibuff + 2, t)

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

    Δt² = Δt*Δt

    return x, y, z, Δt*vx, Δt*vy, Δt*vz, Δt²*ax, Δt²*ay, Δt²*az

end 

"""
    ∂³chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int, Δt, ibuff=1)

Evaluate a sum of Cheybyshev polynomials of the first kind and its three derivatives at `t` 
using a recursive algorithm. It simultenously evalutes the 3 state components. `idx` 
is the index of the starting row (in 0-based notation) in the matrix of coefficients `cₖ`.
`ibuff` is the index of the first free buffer. 

!!! note 
    `x` is a re-work of the actual ascissa value that lies between [-1, 1]

### See Also 
See also [`chebyshev`](@ref), [`∂chebyshev`](@ref) and [`∂²chebyshev`](@ref)

"""
function ∂³chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int, Δt, ibuff=1)

    # Retrieve the work buffer 
    Tₙ    = get_buffer(cache, ibuff, t)
    dTₙ   = get_buffer(cache, ibuff + 1, t)
    ddTₙ  = get_buffer(cache, ibuff + 2, t)
    dddTₙ = get_buffer(cache, ibuff + 3, t)

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

    Δt² = Δt*Δt
    Δt³ = Δt²*Δt

    return x, y, z, Δt*vx, Δt*vy, Δt*vz, Δt²*ax, Δt²*ay, Δt²*az, Δt³*jx, Δt³*jy, Δt³*jz

end 


"""
    ∫chebyshev(cache::InterpCache, cₖ, t::Number, N::Int, Δt, tlen, p₀)

Evaluate the integral of a sum of Cheybyshev polynomials of the first kind using a recursive 
algorithm. It simultenously evalutes the 3 state components. It assumes the Chebyshev polynomials 
up to degree N have already been computed and are stored in the buffer with index `ibuff`. 
`tlen` is the size of the record interval, `Δt` is the timescale factor, and `p₀` is a 
vector containing the position coefficients at the midpoint (i.e., when the integral is 
evaluated at t = 0).

!!! note 
    `x` is a re-work of the actual ascissa value that lies between [-1, 1]
"""
function ∫chebyshev(cache::InterpCache, cₖ, t::Number, N::Int, Δt, tlen, p₀, ibuff)

    # Retrieve the work buffer. The integral buffer is placed in the last position
    Tₙ  = get_buffer(cache, ibuff, t)
    iTₙ = get_buffer(cache, ibuff-1, t)

    # We need to compute this additional Chebyshev polynomial because it is required for 
    # the integrals  
    @inbounds Tₙ[N+1] = 2*t*Tₙ[N] - Tₙ[N-1]

    ix = 1 
    iy = 2 
    iz = 3 

    @inbounds begin 

        iTₙ[1] = t 
        iTₙ[2] = t*t/2 

        x = t*cₖ[ix, 1] + iTₙ[2]*cₖ[ix, 2] 
        y = t*cₖ[iy, 1] + iTₙ[2]*cₖ[iy, 2] 
        z = t*cₖ[iz, 1] + iTₙ[2]*cₖ[iz, 2]

        for j = 4:N+1

            iTₙ[j-1] = (Tₙ[j]/(j-1) - Tₙ[j-2]/(j-3))/2

            # Starting from the integral of T5, the above formula misses a constant of 
            # integration on all the odd polynomials (T5, T7, T9, T11, etc...). The sign of 
            # this constant also alternatively changes between 1 and -1. 
            # Here we are including that constant to ensure that the integral expression is 
            # null when evaluated at t = 0. 

            if isodd(j)
                d = (1/(j-1) + 1/(j-3))/2
                iTₙ[j-1] += ((j-5) % 4 == 0) ? -d : d
            end

            x += iTₙ[j-1]*cₖ[ix, j-1] 
            y += iTₙ[j-1]*cₖ[iy, j-1]
            z += iTₙ[j-1]*cₖ[iz, j-1]

        end

        # Multiplies by tlen/2 because the t coordinate we are using is expressed 
        # as function of the original epoch as t = 2(x - mid)/tlen, so that 
        # dt/dx = 2/tlen 
        x = p₀[1] + x*tlen/Δt/2
        y = p₀[2] + y*tlen/Δt/2
        z = p₀[3] + z*tlen/Δt/2
        
    end

    return x, y, z

end 

"""
    ∫chebyshev(cache::InterpCache, cₖ, t::Number, N::Int, Δt, tlen, p₀)

Evaluate the integral of a sum of Cheybyshev polynomials of the first kind using a recursive 
algorithm. It simultenously evalutes the 3 state components. This function simultaneously 
computes both the Chebyshev polynomials as well as their integrals.
"""
function ∫chebyshev(cache::InterpCache, cₖ, t::Number, N::Int, Δt, tlen, p₀)

    # Retrieve the work buffer. The integral buffer is placed in the last position
    iTₙ = get_buffer(cache, 1, t)
    Tₙ  = get_buffer(cache, 2, t)

    ix = 1 
    iy = 2 
    iz = 3 

    @inbounds begin 
    
        Tₙ[1] = 1 
        Tₙ[2] = t 
        Tₙ[3] = 2*t*t - 1 

        iTₙ[1] = t 
        iTₙ[2] = t*t/2 

        x = t*cₖ[ix, 1] + iTₙ[2]*cₖ[ix, 2] 
        y = t*cₖ[iy, 1] + iTₙ[2]*cₖ[iy, 2] 
        z = t*cₖ[iz, 1] + iTₙ[2]*cₖ[iz, 2]

        for j = 4:N+1

            Tₙ[j] = 2*t*Tₙ[j-1] - Tₙ[j-2]
            iTₙ[j-1] = (Tₙ[j]/(j-1) - Tₙ[j-2]/(j-3))/2

            # Starting from the integral of T5, the above formula misses a constant of 
            # integration on all the odd polynomials (T5, T7, T9, T11, etc...). The sign of 
            # this constant also alternatively changes between 1 and -1. 
            # Here we are including that constant to ensure that the integral expression is 
            # null when evaluated at t = 0. 

            if isodd(j)
                d = (1/(j-1) + 1/(j-3))/2
                iTₙ[j-1] += ((j-5) % 4 == 0) ? -d : d
            end

            x += iTₙ[j-1]*cₖ[ix, j-1] 
            y += iTₙ[j-1]*cₖ[iy, j-1]
            z += iTₙ[j-1]*cₖ[iz, j-1]

        end

        # Multiplies by tlen/2 because the t coordinate we are using is expressed 
        # as function of the original epoch as t = 2(x - mid)/tlen, so that 
        # dt/dx = 2/tlen 
        x = p₀[1] + x*tlen/Δt/2
        y = p₀[2] + y*tlen/Δt/2
        z = p₀[3] + z*tlen/Δt/2
        
    end

    return x, y, z

end 
