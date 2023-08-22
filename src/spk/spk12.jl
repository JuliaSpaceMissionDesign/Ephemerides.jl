
""" 
    SPKSegmentType12(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 12.
"""
function SPKSegmentType12(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader8(daf, desc)
    caches = [SPKSegmentCache8(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType12(header, caches)

end

@inline spk_field(::SPKSegmentType12) = SPK_SEGMENTLIST_MAPPING[12]


function spk_vector3(daf::DAF, seg::SPKSegmentType12, time::Number) 

    index = find_logical_record(header(seg), time)
    get_coefficients!(daf, header(seg), cache(seg), index)

    tbegin = header(seg).tstart + (index-1)*header(seg).tlen 
    Δt = (time - tbegin)/header(seg).tlen + 1

    x = hermite(cache(seg), Δt, 1, header(seg).N, header(seg).tlen)
    y = hermite(cache(seg), Δt, 2, header(seg).N, header(seg).tlen)
    z = hermite(cache(seg), Δt, 3, header(seg).N, header(seg).tlen)

    return SA[x, y, z]

end

function spk_vector6(daf::DAF, seg::SPKSegmentType12, time::Number)

    index = find_logical_record(header(seg), time)
    get_coefficients!(daf, header(seg), cache(seg), index)

    tbegin = header(seg).tstart + (index-1)*header(seg).tlen 
    Δt = (time - tbegin)/header(seg).tlen + 1

    x, vx = ∂hermite(cache(seg), Δt, 1, header(seg).N, header(seg).tlen)
    y, vy = ∂hermite(cache(seg), Δt, 2, header(seg).N, header(seg).tlen)
    z, vz = ∂hermite(cache(seg), Δt, 3, header(seg).N, header(seg).tlen)

    return SA[x, y, z, vx, vy, vz]

end

function spk_vector9(daf::DAF, seg::SPKSegmentType12, time::Number)

    index = find_logical_record(header(seg), time)
    get_coefficients!(daf, header(seg), cache(seg), index)

    tbegin = header(seg).tstart + (index-1)*header(seg).tlen 
    Δt = (time - tbegin)/header(seg).tlen + 1

    x, vx, ax = ∂²hermite(cache(seg), Δt, 1, header(seg).N, header(seg).tlen)
    y, vy, ay = ∂²hermite(cache(seg), Δt, 2, header(seg).N, header(seg).tlen)
    z, vz, az = ∂²hermite(cache(seg), Δt, 3, header(seg).N, header(seg).tlen)

    return SA[x, y, z, vx, vy, vz, ax, ay, az]

end

function spk_vector12(daf::DAF, seg::SPKSegmentType12, time::Number)

    index = find_logical_record(header(seg), time)
    get_coefficients!(daf, header(seg), cache(seg), index)

    tbegin = header(seg).tstart + (index-1)*header(seg).tlen 
    Δt = (time - tbegin)/header(seg).tlen + 1

    x, vx, ax, jx = ∂³hermite(cache(seg), Δt, 1, header(seg).N, header(seg).tlen)
    y, vy, ay, jy = ∂³hermite(cache(seg), Δt, 2, header(seg).N, header(seg).tlen)
    z, vz, az, jz = ∂³hermite(cache(seg), Δt, 3, header(seg).N, header(seg).tlen)

    return SA[x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz]
end


function hermite(cache::SPKSegmentCache8, x::Number, istate::Integer, N::Int, Δt::Number)

    # This function is valid only for equally-spaced polynomials!
    work = get_tmp(cache.work, x)

    # n is the number of points defining the polynomial! 
    # x is a re-work of the ascissa that starts at 1
    # istate is the index of the desired state

    # work is filled in this way [f(x1), df(x1), f(x2), df(x2), ....]
    # We already normalise the derivatives (df/ds = df/dx*dx/ds and dx/ds = h)

    @inbounds for i = 1:N
        work[2i-1] = cache.states[i, istate]
        work[2i]   = cache.states[i, istate+3]*Δt
    end

    # Compute the 2nd column of the interpolation table (N-1 values)
    @inbounds for i = 1:N-1

        c2 = x - i 
        c1 = 1 - c2 

        prev = 2i - 1 
        this = prev + 1 
        next = this + 1

        temp = work[prev] + (x - i)*work[this]
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


function ∂hermite(cache::SPKSegmentCache8, x::Number, istate::Integer, N::Int, Δt::Number)

    # This function is valid only for equally-spaced polynomials
    work  = get_tmp(cache.work, x) 
    dwork = get_tmp(cache.dwork, x)

    # n is the number of points defining the polynomial! 
    # x is a re-work of the ascissa that starts at 1
    # istate is the index of the desired state

    # work is filled in this way [f(x1), df(x1), f(x2), df(x2), ....]
    # We already normalise the derivatives (df/ds = df/dx*dx/ds and dx/ds = h)

    @inbounds for i = 1:N
        work[2i-1] = cache.states[i, istate]
        work[2i]   = cache.states[i, istate+3]*Δt
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

        temp = work[prev] + (x - i)*work[this]
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

function ∂²hermite(cache::SPKSegmentCache8, x::Number, istate::Integer, N::Int, Δt::Number)

    # This function is valid only for equally-spaced polynomials
    work  = get_tmp(cache.work, x) 
    dwork = get_tmp(cache.dwork, x)
    d²work = get_tmp(cache.ddwork, x)
        
    # n is the number of points defining the polynomial! 
    # x is a re-work of the ascissa that starts at 1
    # istate is the index of the desired state

    # work is filled in this way [f(x1), df(x1), f(x2), df(x2), ....]
    # We already normalise the derivatives (df/ds = df/dx*dx/ds and dx/ds = h)

    @inbounds for i = 1:N
        work[2i-1] = cache.states[i, istate]
        work[2i]   = cache.states[i, istate+3]*Δt

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

        temp = work[prev] + (x - i)*work[this]
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

function ∂³hermite(cache::SPKSegmentCache8, x::Number, istate::Integer, N::Int, Δt::Number)

    # This function is valid only for equally-spaced polynomials
    work  = get_tmp(cache.work, x) 
    dwork = get_tmp(cache.dwork, x)
    d²work = get_tmp(cache.ddwork, x)
    d³work = get_tmp(cache.dddwork, x)
        
    # n is the number of points defining the polynomial! 
    # x is a re-work of the ascissa that starts at 1
    # istate is the index of the desired state

    # work is filled in this way [f(x1), df(x1), f(x2), df(x2), ....]
    # We already normalise the derivatives (df/ds = df/dx*dx/ds and dx/ds = h)

    @inbounds for i = 1:N
        work[2i-1] = cache.states[i, istate]
        work[2i]   = cache.states[i, istate+3]*Δt

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

        temp = work[prev] + (x - i)*work[this]
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