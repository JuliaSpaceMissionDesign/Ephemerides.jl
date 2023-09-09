
""" 
    SPKSegmentType13(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 13.
"""
function SPKSegmentType13(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader9(daf, desc)
    caches = [SPKSegmentCache9(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType13(header, caches)

end

@inline spk_field(::SPKSegmentType13) = SPK_SEGMENTLIST_MAPPING[13]


function spk_vector3(daf::DAF, seg::SPKSegmentType13, t::Number) 

    index = find_logical_record(daf, header(seg), t)
    get_coefficients!(daf, header(seg), cache(seg), index)

    x = hermite(cache(seg), t, 1, header(seg).N)
    y = hermite(cache(seg), t, 2, header(seg).N)
    z = hermite(cache(seg), t, 3, header(seg).N)

    return SA[x, y, z]

end

function spk_vector6(daf::DAF, seg::SPKSegmentType13, t::Number)

    index = find_logical_record(daf, header(seg), t)
    get_coefficients!(daf, header(seg), cache(seg), index)

    x, vx = ∂hermite(cache(seg), t, 1, header(seg).N)
    y, vy = ∂hermite(cache(seg), t, 2, header(seg).N)
    z, vz = ∂hermite(cache(seg), t, 3, header(seg).N)

    return SA[x, y, z, vx, vy, vz]

end

function spk_vector9(daf::DAF, seg::SPKSegmentType13, t::Number)

    index = find_logical_record(daf, header(seg), t)
    get_coefficients!(daf, header(seg), cache(seg), index)

    x, vx, ax = ∂²hermite(cache(seg), t, 1, header(seg).N)
    y, vy, ay = ∂²hermite(cache(seg), t, 2, header(seg).N)
    z, vz, az = ∂²hermite(cache(seg), t, 3, header(seg).N)

    return SA[x, y, z, vx, vy, vz, ax, ay, az]

end

function spk_vector12(daf::DAF, seg::SPKSegmentType13, t::Number)

    index = find_logical_record(daf, header(seg), t)
    get_coefficients!(daf, header(seg), cache(seg), index)

    x, vx, ax, jx = ∂³hermite(cache(seg), t, 1, header(seg).N)
    y, vy, ay, jy = ∂³hermite(cache(seg), t, 2, header(seg).N)
    z, vz, az, jz = ∂³hermite(cache(seg), t, 3, header(seg).N)

    return SA[x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz]
end


function hermite(cache::SPKSegmentCache9, x::Number, istate::Integer, N::Int)

    # This function is also valid for unequal time steps 
    # istate is the index of the desired state
    
    work = get_tmp(cache.work, x)

    # work is filled in this way [f(x1), df(x1), f(x2), df(x2), ....]
    @inbounds for i = 1:N
        work[2i-1] = cache.states[i, istate]
        work[2i]   = cache.states[i, istate+3]
    end

    # Compute the 2nd column of the interpolation table (N-1 values)
    @inbounds for i = 1:N-1

        c1 = cache.epochs[i+1] - x 
        c2 = x - cache.epochs[i]
        d = c1 + c2

        prev = 2i - 1 
        this = prev + 1 
        next = this + 1

        temp = work[prev] + c2*work[this]
        work[this] = (c1*work[prev] + c2*work[next])/d
        work[prev] = temp
        
    end

    @inbounds work[2N-1] = work[2N-1] + work[2N]*(x-cache.epochs[N])

    # Compute columns 3 through 2N of the table 
    @inbounds for j = 2:2N-1
        for i = 1:2N-j

            xi = (i+1) ÷ 2 
            xij = (i+j+1) ÷ 2 

            c1 = cache.epochs[xij] - x 
            c2 = x - cache.epochs[xi]
            d = c1 + c2

            work[i] = (c1*work[i] + c2*work[i+1])/d
        end
    end

    @inbounds return work[1]

end 


function ∂hermite(cache::SPKSegmentCache9, x::Number, istate::Integer, N::Int)

    # This function is also valid for unequal time steps 
    # istate is the index of the desired state

    work  = get_tmp(cache.work, x) 
    dwork = get_tmp(cache.dwork, x)

    # work is filled in this way [f(x1), df(x1), f(x2), df(x2), ....]
    @inbounds for i = 1:N
        work[2i-1] = cache.states[i, istate]
        work[2i]   = cache.states[i, istate+3]
    end

    # Compute the 2nd column of the interpolation table (N-1 values)
    @inbounds for i = 1:N-1

        c1 = cache.epochs[i+1] - x 
        c2 = x - cache.epochs[i]
        d = c1 + c2

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
        work[2N-1] = work[2N-1] + work[2N]*(x-cache.epochs[N])
        dwork[2N-1] = work[2N]
    end

    # Compute columns 3 through 2N of the table 
    @inbounds for j = 2:2N-1
        for i = 1:2N-j

            xi = (i+1) ÷ 2 
            xij = (i+j+1) ÷ 2 

            c1 = cache.epochs[xij] - x 
            c2 = x - cache.epochs[xi] 
            d = c1 + c2

            dwork[i] = (c1*dwork[i] + c2*dwork[i+1] + work[i+1] - work[i])/d
            work[i] = (c1*work[i] + c2*work[i+1])/d
        end
    end

    @inbounds return work[1], dwork[1]

end 

function ∂²hermite(cache::SPKSegmentCache9, x::Number, istate::Integer, N::Int)

    # This function is also valid for unequal time steps 
    # istate is the index of the desired state

    work  = get_tmp(cache.work, x) 
    dwork = get_tmp(cache.dwork, x)
    d²work = get_tmp(cache.ddwork, x)
        
    # work is filled in this way [f(x1), df(x1), f(x2), df(x2), ....]
    @inbounds for i = 1:N
        work[2i-1] = cache.states[i, istate]
        work[2i]   = cache.states[i, istate+3]

        d²work[2i-1] = 0
        d²work[2i] = 0
    end

    # Compute the 2nd column of the interpolation table (N-1 values)
    @inbounds for i = 1:N-1

        c1 = cache.epochs[i+1] - x 
        c2 = x - cache.epochs[i]
        d = c1 + c2

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
        work[2N-1] = work[2N-1] + work[2N]*(x-cache.epochs[N])
        dwork[2N-1] = work[2N]
    end

    # Compute columns 3 through 2N of the table 
    @inbounds for j = 2:2N-1
        for i = 1:2N-j

            xi = (i+1) ÷ 2 
            xij = (i+j+1) ÷ 2 

            c1 = cache.epochs[xij] - x 
            c2 = x - cache.epochs[xi]
            d = c1 + c2 

            d²work[i] = (c1*d²work[i] + c2*d²work[i+1] + 2*(dwork[i+1] - dwork[i]))/d
            dwork[i]  = (c1*dwork[i]  + c2*dwork[i+1] + work[i+1] - work[i])/d
            work[i]   = (c1*work[i]   + c2*work[i+1])/d
        end
    end

    @inbounds return work[1], dwork[1], d²work[1]

end 

function ∂³hermite(cache::SPKSegmentCache9, x::Number, istate::Integer, N::Int)

    # This function is also valid for unequal time steps 
    # istate is the index of the desired state

    work  = get_tmp(cache.work, x) 
    dwork = get_tmp(cache.dwork, x)
    d²work = get_tmp(cache.ddwork, x)
    d³work = get_tmp(cache.dddwork, x)
        
    # work is filled in this way [f(x1), df(x1), f(x2), df(x2), ....]
    @inbounds for i = 1:N
        work[2i-1] = cache.states[i, istate]
        work[2i]   = cache.states[i, istate+3]

        d²work[2i-1] = 0
        d²work[2i] = 0

        d³work[2i-1] = 0
        d³work[2i] = 0
    end

    # Compute the 2nd column of the interpolation table (N-1 values)
    @inbounds for i = 1:N-1

        c1 = cache.epochs[i+1] - x 
        c2 = x - cache.epochs[i]
        d = c1 + c2

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
        work[2N-1] = work[2N-1] + work[2N]*(x-cache.epochs[N])
        dwork[2N-1] = work[2N]
    end

    # Compute columns 3 through 2N of the table 
    @inbounds for j = 2:2N-1
        for i = 1:2N-j

            xi = (i+1) ÷ 2 
            xij = (i+j+1) ÷ 2 

            c1 = cache.epochs[xij] - x 
            c2 = x - cache.epochs[xi] 
            d = c1 + c2
            
            d³work[i] = (c1*d³work[i] + c2*d³work[i+1] + 3*(d²work[i+1] - d²work[i]))/d
            d²work[i] = (c1*d²work[i] + c2*d²work[i+1] + 2*(dwork[i+1] - dwork[i]))/d
            dwork[i]  = (c1*dwork[i]  + c2*dwork[i+1] + work[i+1] - work[i])/d
            work[i]   = (c1*work[i]   + c2*work[i+1])/d
        end
    end

    @inbounds return work[1], dwork[1], d²work[1],  d³work[1]

end 