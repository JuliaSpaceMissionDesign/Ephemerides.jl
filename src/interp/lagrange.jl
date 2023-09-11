
"""
    lagrange(cache::InterpCache, states, x, idx::Int, N::Int)

Recursively evaluate a Lagrange polynomial at `x` by using Neville's algorithm. This 
function is valid only for equally-spaced polynomials. `idx` is the index of the desired 
state and `N` is the number of coefficients of the polynomial. 

!!! note 
    `x` is a re-work of the actual ascissa value that starts at 1 

### See Also 
See also [`∂lagrange`](@ref) and [`∂²lagrange`](@ref).
"""
function lagrange(cache::InterpCache, states, x::Number, idx::Int, N::Int)

    # Retrieve the work buffer 
    work = get_buffer(cache, 1, x)

    # Initialise the buffer
    @inbounds for i = 1:N 
        work[i] = states[i, idx]
    end

    # Compute the polynomials using a recursive relationship
    @inbounds for j = 1:N-1 
        for i = 1:N-j 
            c1 = i + j - x
            c2 = x - i 

            work[i] = (c1*work[i] + c2*work[i+1])/j
        end
    end

    return work[1]

end 

"""
    lagrange(cache::InterpCache, states, epochs, x, idx::Int, N::Int)

Recursively evaluate a Lagrange polynomial at `x` by using Neville's algorithm. This 
function handles unequally-spaced polynomials, where the coefficients in `states` are 
interpolated at `epochs`.
"""
function lagrange(cache::InterpCache, states, epochs, x::Number, idx::Int, N::Int)

    # Retrieve the work buffer 
    work = get_buffer(cache, 1, x)
    
    # Initialise the buffer
    @inbounds for i = 1:N 
        work[i] = states[i, idx]
    end

    # Compute the polynomials using a recursive relationship
    @inbounds for j = 1:N-1 
        for i = 1:N-j 

            c1 = x - epochs[i+j]
            c2 = epochs[i] - x
            d = epochs[i] - epochs[i+j]

            work[i] = (c1*work[i] + c2*work[i+1])/d
        end
    end

    return work[1]

end 


"""
    ∂lagrange(cache::InterpCache, states, x, idx::Int, N::Int, Δt::Number)

Recursively evaluate a Lagrange polynomial and its derivative at `x` by using 
Neville's algorithm. This function is valid only for equally-spaced polynomials. `idx` 
is the index of the desired state, `N` is the number of coefficients of the polynomial and 
`Δt`is the length of the polynomial interval

!!! note 
    `x` is a re-work of the actual ascissa value that starts at 1

### See Also 
See also [`lagrange`](@ref) and [`∂²lagrange`](@ref)
"""
function ∂lagrange(cache::InterpCache, states, x::Number, idx::Int, N::Int, Δt::Number)

    # Retrieve the work buffers
    work  = get_buffer(cache, 1, x)
    dwork = get_buffer(cache, 2, x)

    # Initialise the buffers
    @inbounds for i = 1:N 
        work[i] = states[i, idx]
        dwork[i] = 0.0
    end

    # Precompute abscissa derivatives
    dc2 = 1/Δt
    dc1 = -dc2

    # Compute the polynomials using a recursive relationship
    @inbounds for j = 1:N-1 
        for i = 1:N-j 
            c1 = i + j - x
            c2 = x - i 

            dwork[i] = (dc1*work[i]   + c1*dwork[i] + dc2*work[i+1] + c2*dwork[i+1])/j
            work[i] = (c1*work[i] + c2*work[i+1])/j

        end
    end

    return work[1], dwork[1]

end 

"""
    ∂lagrange(cache::InterpCache, states, epochs, x, idx::Int, N::Int)

Recursively evaluate a Lagrange polynomial and its derivative at `x` by using 
Neville's algorithm. This function handles unequally-spaced polynomials, where the 
coefficients in `states` are interpolated at `epochs`.
"""
function ∂lagrange(cache::InterpCache, states, epochs, x::Number, idx::Int, N::Int)

    # Retrieve the work buffers
    work  = get_buffer(cache, 1, x)
    dwork = get_buffer(cache, 2, x)

    # Initialise the buffers
    @inbounds for i = 1:N 
        work[i] = states[i, idx]
        dwork[i] = 0.0
    end

    # Compute the polynomials using a recursive relationship
    @inbounds for j = 1:N-1 
        for i = 1:N-j 

            c1 = x - epochs[i+j]
            c2 = epochs[i] - x            
            d = epochs[i] - epochs[i+j]

            dwork[i] = (work[i] + c1*dwork[i] - work[i+1] + c2*dwork[i+1])/d
            work[i] = (c1*work[i] + c2*work[i+1])/d

        end
    end

    return work[1], dwork[1]

end 


"""
    ∂²lagrange(cache::InterpCache, states, x, idx::Int, N::Int, Δt::Number)

Recursively evaluate a Lagrange polynomial and its two derivatives at `x` by using 
Neville's algorithm. This function is valid only for equally-spaced polynomials. `idx` 
is the index of the desired state, `N` is the number of coefficients of the polynomial and 
`Δt`is the length of the polynomial interval

!!! note 
    `x` is a re-work of the actual ascissa value that starts at 1

### See Also 
See also [`lagrange`](@ref) and [`∂lagrange`](@ref)
"""
function ∂²lagrange(cache::InterpCache, states, x::Number, idx::Int, N::Int, Δt::Number)

    # Retrieve the work buffers
    work   = get_buffer(cache, 1, x)
    dwork  = get_buffer(cache, 2, x)
    ddwork = get_buffer(cache, 3, x)

    # Initialise the buffers
    @inbounds for i = 1:N 
        work[i] = states[i, idx]
        dwork[i] = 0.0
        ddwork[i] = 0.0
    end

    # Precompute abscissa derivatives
    dc2 = 1/Δt
    dc1 = -dc2

    # Compute the polynomials using a recursive relationship
    @inbounds for j = 1:N-1 
        for i = 1:N-j 
            c1 = i + j - x
            c2 = x - i 

            ddwork[i] = (2*dc1*dwork[i]   + c1*ddwork[i] + 
                         2*dc2*dwork[i+1] + c2*ddwork[i+1])/j

            dwork[i] = (dc1*work[i]   + c1*dwork[i] + dc2*work[i+1] + c2*dwork[i+1])/j
            work[i] = (c1*work[i] + c2*work[i+1])/j

        end
    end

    return work[1], dwork[1], ddwork[1]

end 

"""
    ∂²lagrange(cache::InterpCache, states, epochs, x, idx::Int, N::Int)

Recursively evaluate a Lagrange polynomial and its two derivatives at `x` by using 
Neville's algorithm. This function handles unequally-spaced polynomials, where the 
coefficients in `states` are interpolated at `epochs`.

### See Also 
See also [`lagrange`](@ref) and [`∂lagrange`](@ref)
"""
function ∂²lagrange(cache::InterpCache, states, epochs, x::Number, idx::Int, N::Int)

    # Retrieve the work buffers
    work   = get_buffer(cache, 1, x)
    dwork  = get_buffer(cache, 2, x)
    ddwork = get_buffer(cache, 3, x)

    # Initialise the buffers
    @inbounds for i = 1:N 
        work[i] = states[i, idx]
        dwork[i] = 0.0
        ddwork[i] = 0.0
    end

    # Compute the polynomials using a recursive relationship
    @inbounds for j = 1:N-1 
        for i = 1:N-j 

            c1 = x - epochs[i+j]
            c2 = epochs[i] - x            
            d = epochs[i] - epochs[i+j]

            ddwork[i] = (2*dwork[i] + c1*ddwork[i] - 2*dwork[i+1] + c2*ddwork[i+1])/d
            dwork[i] = (work[i] + c1*dwork[i] - work[i+1] + c2*dwork[i+1])/d
            work[i] = (c1*work[i] + c2*work[i+1])/d

        end
    end

    return work[1], dwork[1], ddwork[1]

end 
