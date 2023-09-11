
"""
    hermite(cache::InterpCache, states, x, idx::Int, N::Int, Δt::Number)

Evalute a Hermite polynomial at `x` using a recursive algorithm. This function is valid 
only for equally-spaced polynomials. `idx` is the index of the desired state, `N` 
is the number of coefficients of the polynomial and `Δt` is the length of the polynomial 
interval.

!!! note 
    `x` is a re-work of the actual ascissa value that starts at 1 

### See Also 
See also [`∂hermite`](@ref), [`∂²hermite`](@ref) and [`∂³hermite`](@ref).
"""
function hermite(cache::InterpCache, states, x::Number, idx::Int, N::Int, Δt::Number)

    # Retrieve the work buffer 
    work = get_buffer(cache, 1, x)

    # Initialise the buffer in this way [f(x1), df(x1), f(x2), df(x2), ....]
    # We already normalise the derivatives (df/ds = df/dx*dx/ds and dx/ds = h)
    @inbounds for i = 1:N
        work[2i-1] = states[i, idx]
        work[2i]   = states[i, idx+3]*Δt
    end

    # Compute the 2nd column of the interpolation table (N-1 values)
    @inbounds for i = 1:N-1

        c2 = x - i 
        c1 = 1 - c2 

        prev = 2i - 1 
        this = prev + 1 
        next = this + 1

        temp = work[prev] + c2*work[this]
        work[this] = c1*work[prev] + c2*work[next]
        work[prev] = temp
        
    end

    @inbounds work[2N-1] = work[2N-1] + work[2N]*(x-N)

    # Compute columns 3 through 2N of the table 
    @inbounds for j = 2:2N-1
        for i = 1:2N-j

            xi = (i+1) ÷ 2 
            xij = (i+j+1) ÷ 2 

            c1 = xij - x 
            c2 = x - xi 
            
            denom = xij - xi 

            work[i] = (c1*work[i] + c2*work[i+1])/denom
        end
    end

    @inbounds return work[1]

end 

"""
    hermite(cache::InterpCache, states, epochs, x, idx::Int, N::Int)

Evalute a Hermite polynomial at `x` using a recursive algorithm. This function handles 
unequally-spaced polynomials, where the coefficients in `states` are interpolated 
at `epochs`.
"""
function hermite(cache::InterpCache, states, epochs, x::Number, idx::Int, N::Int)

    # Retrieve the work buffer 
    work = get_buffer(cache, 1, x)

    # Initialise the buffer in this way [f(x1), df(x1), f(x2), df(x2), ....]
    @inbounds for i = 1:N
        work[2i-1] = states[i, idx]
        work[2i]   = states[i, idx+3]
    end

    # Compute the 2nd column of the interpolation table (N-1 values)
    @inbounds for i = 1:N-1

        c1 = epochs[i+1] - x 
        c2 = x - epochs[i]
        d = epochs[i+1] - epochs[i]

        prev = 2i - 1 
        this = prev + 1 
        next = this + 1

        temp = work[prev] + c2*work[this]
        work[this] = (c1*work[prev] + c2*work[next])/d
        work[prev] = temp
        
    end

    @inbounds work[2N-1] = work[2N-1] + work[2N]*(x-epochs[N])

    # Compute columns 3 through 2N of the table 
    @inbounds for j = 2:2N-1
        for i = 1:2N-j

            xi = (i+1) ÷ 2 
            xij = (i+j+1) ÷ 2 

            c1 = epochs[xij] - x 
            c2 = x - epochs[xi]
            d = epochs[xij] - epochs[xi]

            work[i] = (c1*work[i] + c2*work[i+1])/d
        end
    end

    @inbounds return work[1]

end 


"""
    ∂hermite(cache::InterpCache, states, x, idx::Int, N::Int, Δt::Number)

Evalute a Hermite polynomial and its derivative at `x` using a recursive algorithm. 
This function is valid only for equally-spaced polynomials. `idx` is the index of the 
desired state, `N` is the number of coefficients of the polynomial and `Δt` is the length 
of the polynomial interval.

!!! note 
    `x` is a re-work of the actual ascissa value that starts at 1 

### See Also 
See also [`hermite`](@ref), [`∂²hermite`](@ref) and [`∂³hermite`](@ref).
"""
function ∂hermite(cache::InterpCache, states, x::Number, idx::Int, N::Int, Δt::Number)

    # Retrieve the work buffers
    work  = get_buffer(cache, 1, x)
    dwork = get_buffer(cache, 2, x)

    # Initialise the buffers in this way [f(x1), df(x1), f(x2), df(x2), ....]
    # We already normalise the derivatives (df/ds = df/dx*dx/ds and dx/ds = h)
    @inbounds for i = 1:N
        work[2i-1] = states[i, idx]
        work[2i]   = states[i, idx+3]*Δt
    end

    # Compute the 2nd column of the interpolation table (N-1 values)
    @inbounds for i = 1:N-1

        c2 = x - i 
        c1 = 1 - c2 

        prev = 2i - 1 
        this = prev + 1 
        next = this + 1

        # Compute the derivative first 
        dwork[prev] = work[this]
        dwork[this] = work[next] - work[prev]

        temp = work[prev] + c2*work[this]
        work[this] = c1*work[prev] + c2*work[next]
        work[prev] = temp
        
    end

    @inbounds begin 
        work[2N-1] = work[2N-1] + work[2N]*(x-N)
        dwork[2N-1] = work[2N]
    end

    # Compute columns 3 through 2N of the table 
    @inbounds for j = 2:2N-1
        for i = 1:2N-j

            xi = (i+1) ÷ 2 
            xij = (i+j+1) ÷ 2 

            c1 = xij - x 
            c2 = x - xi 
            
            denom = xij - xi 

            dwork[i] = (c1*dwork[i] + c2*dwork[i+1] + work[i+1] - work[i])/denom
            work[i] = (c1*work[i] + c2*work[i+1])/denom
        end
    end

    @inbounds return work[1], dwork[1]/Δt

end 

"""
    ∂hermite(cache::InterpCache, states, epochs, x, idx::Int, N::Int)

Evalute a Hermite polynomial and its derivative at `x` using a recursive algorithm. 
This function handles unequally-spaced polynomials, where the coefficients in `states` 
are interpolated at `epochs`.
"""
function ∂hermite(cache::InterpCache, states, epochs, x::Number, idx::Int, N::Int)

    # Retrieve the work buffers
    work  = get_buffer(cache, 1, x)
    dwork = get_buffer(cache, 2, x)

    # Initialise the buffers in this way [f(x1), df(x1), f(x2), df(x2), ....]
    @inbounds for i = 1:N
        work[2i-1] = states[i, idx]
        work[2i]   = states[i, idx+3]
    end

    # Compute the 2nd column of the interpolation table (N-1 values)
    @inbounds for i = 1:N-1

        c1 = epochs[i+1] - x 
        c2 = x - epochs[i]
        d = epochs[i+1] - epochs[i]

        prev = 2i - 1 
        this = prev + 1 
        next = this + 1

        # Compute the derivative first 
        dwork[prev] = work[this]
        dwork[this] = (work[next] - work[prev])/d

        temp = work[prev] + c2*work[this]
        work[this] = (c1*work[prev] + c2*work[next])/d
        work[prev] = temp
        
    end

    @inbounds begin 
        work[2N-1] = work[2N-1] + work[2N]*(x-epochs[N])
        dwork[2N-1] = work[2N]
    end

    # Compute columns 3 through 2N of the table 
    @inbounds for j = 2:2N-1
        for i = 1:2N-j

            xi = (i+1) ÷ 2 
            xij = (i+j+1) ÷ 2 

            c1 = epochs[xij] - x 
            c2 = x - epochs[xi] 
            d = epochs[xij] - epochs[xi]

            dwork[i] = (c1*dwork[i] + c2*dwork[i+1] + work[i+1] - work[i])/d
            work[i] = (c1*work[i] + c2*work[i+1])/d
        end
    end

    @inbounds return work[1], dwork[1]

end 


"""
    ∂²hermite(cache::InterpCache, states, x, idx::Int, N::Int, Δt::Number)

Evalute a Hermite polynomial and its two derivatives at `x` using a recursive algorithm. 
This function is valid only for equally-spaced polynomials. `idx` is the index of the 
desired state, `N` is the number of coefficients of the polynomial and `Δt` is the length 
of the polynomial interval.

!!! note 
    `x` is a re-work of the actual ascissa value that starts at 1 

### See Also 
See also [`hermite`](@ref), [`∂hermite`](@ref) and [`∂³hermite`](@ref).
"""
function ∂²hermite(cache::InterpCache, states, x::Number, idx::Int, N::Int, Δt::Number)

    # Retrieve the work buffers
    work   = get_buffer(cache, 1, x)
    dwork  = get_buffer(cache, 2, x)
    d²work = get_buffer(cache, 3, x)
        
    # Initialise the buffers in this way [f(x1), df(x1), f(x2), df(x2), ....]
    # We already normalise the derivatives (df/ds = df/dx*dx/ds and dx/ds = h)
    @inbounds for i = 1:N
        work[2i-1] = states[i, idx]
        work[2i]   = states[i, idx+3]*Δt

        d²work[2i-1] = 0
        d²work[2i] = 0
    end

    # Compute the 2nd column of the interpolation table (N-1 values)
    @inbounds for i = 1:N-1

        c2 = x - i 
        c1 = 1 - c2 

        prev = 2i - 1 
        this = prev + 1 
        next = this + 1

        # Compute the derivatives first 
        dwork[prev] = work[this]
        dwork[this] = work[next] - work[prev]

        temp = work[prev] + c2*work[this]
        work[this] = c1*work[prev] + c2*work[next]
        work[prev] = temp
        
    end

    @inbounds begin 
        work[2N-1] = work[2N-1] + work[2N]*(x-N)
        dwork[2N-1] = work[2N]
    end

    # Compute columns 3 through 2N of the table 
    @inbounds for j = 2:2N-1
        for i = 1:2N-j

            xi = (i+1) ÷ 2 
            xij = (i+j+1) ÷ 2 

            c1 = xij - x 
            c2 = x - xi 
            
            denom = xij - xi 

            d²work[i] = (c1*d²work[i] + c2*d²work[i+1] + 2*(dwork[i+1] - dwork[i]))/denom
            dwork[i]  = (c1*dwork[i]  + c2*dwork[i+1] + work[i+1] - work[i])/denom
            work[i]   = (c1*work[i]   + c2*work[i+1])/denom
        end
    end

    @inbounds return work[1], dwork[1]/Δt, d²work[1]/Δt^2

end 

"""
    ∂²hermite(cache::InterpCache, states, epochs, x, idx::Int, N::Int)

Evalute a Hermite polynomial and its two derivatives at `x` using a recursive algorithm. 
This function handles unequally-spaced polynomials, where the coefficients in `states` 
are interpolated at `epochs`.
"""
function ∂²hermite(cache::InterpCache, states, epochs, x::Number, idx::Int, N::Int)

    # Retrieve the work buffers
    work   = get_buffer(cache, 1, x)
    dwork  = get_buffer(cache, 2, x)
    d²work = get_buffer(cache, 3, x)
        
    # Initialise the buffers in this way [f(x1), df(x1), f(x2), df(x2), ....]
    @inbounds for i = 1:N
        work[2i-1] = states[i, idx]
        work[2i]   = states[i, idx+3]

        d²work[2i-1] = 0
        d²work[2i] = 0
    end

    # Compute the 2nd column of the interpolation table (N-1 values)
    @inbounds for i = 1:N-1

        c1 = epochs[i+1] - x 
        c2 = x - epochs[i]
        d = epochs[i+1] - epochs[i]

        prev = 2i - 1 
        this = prev + 1 
        next = this + 1

        # Compute the derivatives first 
        dwork[prev] = work[this]
        dwork[this] = (work[next] - work[prev])/d

        temp = work[prev] + c2*work[this]
        work[this] = (c1*work[prev] + c2*work[next])/d
        work[prev] = temp
        
    end

    @inbounds begin 
        work[2N-1] = work[2N-1] + work[2N]*(x-epochs[N])
        dwork[2N-1] = work[2N]
    end

    # Compute columns 3 through 2N of the table 
    @inbounds for j = 2:2N-1
        for i = 1:2N-j

            xi = (i+1) ÷ 2 
            xij = (i+j+1) ÷ 2 

            c1 = epochs[xij] - x 
            c2 = x - epochs[xi]
            d = epochs[xij] - epochs[xi]

            d²work[i] = (c1*d²work[i] + c2*d²work[i+1] + 2*(dwork[i+1] - dwork[i]))/d
            dwork[i]  = (c1*dwork[i]  + c2*dwork[i+1] + work[i+1] - work[i])/d
            work[i]   = (c1*work[i]   + c2*work[i+1])/d
        end
    end

    @inbounds return work[1], dwork[1], d²work[1]

end 


"""
    ∂³hermite(cache::InterpCache, states, x, idx::Int, N::Int, Δt::Number)

Evalute a Hermite polynomial and its three derivatives at `x` using a recursive algorithm. 
This function is valid only for equally-spaced polynomials. `idx` is the index of the 
desired state, `N` is the number of coefficients of the polynomial and `Δt` is the length 
of the polynomial interval.

!!! note 
    `x` is a re-work of the actual ascissa value that starts at 1 

### See Also 
See also [`hermite`](@ref), [`∂hermite`](@ref) and [`∂²hermite`](@ref).
"""
function ∂³hermite(cache::InterpCache, states, x::Number, idx::Int, N::Int, Δt::Number)

    # Retrieve the work buffers
    work   = get_buffer(cache, 1, x)
    dwork  = get_buffer(cache, 2, x)
    d²work = get_buffer(cache, 3, x)
    d³work = get_buffer(cache, 4, x)
        
    # Initialise the buffers in this way [f(x1), df(x1), f(x2), df(x2), ....]
    # We already normalise the derivatives (df/ds = df/dx*dx/ds and dx/ds = h)
    @inbounds for i = 1:N
        work[2i-1] = states[i, idx]
        work[2i]   = states[i, idx+3]*Δt

        d²work[2i-1] = 0
        d²work[2i] = 0

        d³work[2i-1] = 0
        d³work[2i] = 0
    end

    # Compute the 2nd column of the interpolation table (N-1 values)
    @inbounds for i = 1:N-1

        c2 = x - i 
        c1 = 1 - c2 

        prev = 2i - 1 
        this = prev + 1 
        next = this + 1

        # Compute the derivatives first 
        dwork[prev] = work[this]
        dwork[this] = work[next] - work[prev]

        temp = work[prev] + c2*work[this]
        work[this] = c1*work[prev] + c2*work[next]
        work[prev] = temp
        
    end

    @inbounds begin 
        work[2N-1] = work[2N-1] + work[2N]*(x-N)
        dwork[2N-1] = work[2N]
    end

    # Compute columns 3 through 2N of the table 
    @inbounds for j = 2:2N-1
        for i = 1:2N-j

            xi = (i+1) ÷ 2 
            xij = (i+j+1) ÷ 2 

            c1 = xij - x 
            c2 = x - xi 
            
            denom = xij - xi 

            d³work[i] = (c1*d³work[i] + c2*d³work[i+1] + 3*(d²work[i+1] - d²work[i]))/denom
            d²work[i] = (c1*d²work[i] + c2*d²work[i+1] + 2*(dwork[i+1] - dwork[i]))/denom
            dwork[i]  = (c1*dwork[i]  + c2*dwork[i+1] + work[i+1] - work[i])/denom
            work[i]   = (c1*work[i]   + c2*work[i+1])/denom
        end
    end

    Δt² = Δt*Δt
    Δt³ = Δt²*Δt

    @inbounds return work[1], dwork[1]/Δt, d²work[1]/Δt²,  d³work[1]/Δt³

end 

"""
    ∂³hermite(cache::InterpCache, states, epochs, x, idx::Int, N::Int)

Evalute a Hermite polynomial and its three derivatives at `x` using a recursive algorithm. 
This function handles unequally-spaced polynomials, where the coefficients in `states` 
are interpolated at `epochs`.
"""
function ∂³hermite(cache::InterpCache, states, epochs, x::Number, idx::Int, N::Int)

    # Retrieve the work buffers
    work   = get_buffer(cache, 1, x)
    dwork  = get_buffer(cache, 2, x)
    d²work = get_buffer(cache, 3, x)
    d³work = get_buffer(cache, 4, x)
        
    # Initialise the buffers in this way [f(x1), df(x1), f(x2), df(x2), ....]
    @inbounds for i = 1:N
        work[2i-1] = states[i, idx]
        work[2i]   = states[i, idx+3]

        d²work[2i-1] = 0
        d²work[2i] = 0

        d³work[2i-1] = 0
        d³work[2i] = 0
    end

    # Compute the 2nd column of the interpolation table (N-1 values)
    @inbounds for i = 1:N-1

        c1 = epochs[i+1] - x 
        c2 = x - epochs[i]
        d = epochs[i+1] - epochs[i]

        prev = 2i - 1 
        this = prev + 1 
        next = this + 1

        # Compute the derivatives first 
        dwork[prev] = work[this]
        dwork[this] = (work[next] - work[prev])/d

        temp = work[prev] + c2*work[this]
        work[this] = (c1*work[prev] + c2*work[next])/d
        work[prev] = temp
        
    end

    @inbounds begin 
        work[2N-1] = work[2N-1] + work[2N]*(x-epochs[N])
        dwork[2N-1] = work[2N]
    end

    # Compute columns 3 through 2N of the table 
    @inbounds for j = 2:2N-1
        for i = 1:2N-j

            xi = (i+1) ÷ 2 
            xij = (i+j+1) ÷ 2 

            c1 = epochs[xij] - x 
            c2 = x - epochs[xi] 
            d = epochs[xij] - epochs[xi]
            
            d³work[i] = (c1*d³work[i] + c2*d³work[i+1] + 3*(d²work[i+1] - d²work[i]))/d
            d²work[i] = (c1*d²work[i] + c2*d²work[i+1] + 2*(dwork[i+1] - dwork[i]))/d
            dwork[i]  = (c1*dwork[i]  + c2*dwork[i+1] + work[i+1] - work[i])/d
            work[i]   = (c1*work[i]   + c2*work[i+1])/d
        end
    end

    @inbounds return work[1], dwork[1], d²work[1],  d³work[1]

end 