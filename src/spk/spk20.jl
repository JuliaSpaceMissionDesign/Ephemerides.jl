
""" 
    SPKSegmentHeader20(daf::DAF, desc::DAFSegmentDescriptor)

Create the segment header for an SPK segment of type 20
"""
function SPKSegmentHeader20(daf::DAF, desc::DAFSegmentDescriptor)

    i0 = 8*(final_address(desc)-7)

    # Length and time scales in km and seconds
    dscale = get_float(array(daf), i0, endian(daf))
    tscale = get_float(array(daf), i0 + 8, endian(daf))

    # Everything is now transformed in seconds past J2000 and kms 
    DJ2000, D2S = 2451545, 86400

    # Integer and fractional part of the initial juliad date
    initjd = get_float(array(daf), i0 + 16, endian(daf))
    initfr = get_float(array(daf), i0 + 24, endian(daf))

    # Start epoch and interval length (in seconds)
    tstart = ((initjd - DJ2000) + initfr)*D2S
    tlen = get_float(array(daf), i0 + 32, endian(daf))*D2S

    # Number of elements in each record
    rsize = Int(get_float(array(daf), i0 + 40, endian(daf)))

    # Byte size of each logical record
    recsize = 8*rsize 

    # Polynomial degree 
    order = (rsize - 3) ÷ 3 - 1

    # Polynomial group size (number of coefficients required for the interpolation)
    N = order + 1

    # Number of records 
    n = Int(get_float(array(daf), i0 + 48, endian(daf)))

    # Initial segment address 
    iaa = initial_address(desc)

    SPKSegmentHeader20(dscale, tscale, tstart, tlen, recsize, order, N, n, iaa)
end

""" 
    SPKSegmentCache20(spkhead::SPKSegmentHeader2)

Initialise the cache for an SPK segment of type 20.
"""
function SPKSegmentCache20(head::SPKSegmentHeader20) 
    SPKSegmentCache20(
        MVector(-1), 
        MVector(0.0, 0.0, 0.0), 
        zeros(3, max(3, head.N)),
        InterpCache{Float64}(4, max(3, head.N+1))
        )
end

""" 
    SPKSegmentType20(daf::DAF, desc::DAFSegmentDescriptor)

Create the object representing an SPK segment of type 2.
"""
function SPKSegmentType20(daf::DAF, desc::DAFSegmentDescriptor)

    # Initialise the segment header and cache
    header = SPKSegmentHeader20(daf, desc)
    caches = [SPKSegmentCache20(header) for _ in 1:Threads.nthreads()]

    SPKSegmentType20(header, caches)

end

@inline spk_field(::SPKSegmentType20) = SPK_SEGMENTLIST_MAPPING[20]

function spk_vector3(daf::DAF, seg::SPKSegmentType20, time::Number) 

    head = header(seg)
    data = cache(seg)


end


function spk_vector6(daf::DAF, seg::SPKSegmentType20, time::Number)

    head = header(seg)
    data = cache(seg)

    # Find the logical record containing the Chebyshev coefficients at `time`
    index = find_logical_record(head, time)
    get_coefficients!(daf, head, data, index)

    # Normalise the time argument between [-1, 1]
    t = chebyshev_time(head, time, index)

    _ = ∫chebyshev(
        data.buff, data.A, t, 0, head.N, head.dscale, head.tscale, head.tlen, data.p
    )

    x, y, z, vx, vy, vz = ∫cheb(
        data.buff, data.A, t, head.N, head.dscale, head.tscale, head.tlen, data.p
    )

    return SA[x, y, z, vx, vy, vz]

end

function spk_vector9(daf::DAF, seg::SPKSegmentType20, time::Number)

    head = header(seg)
    data = cache(seg)

    # TODO: implement me

end

function spk_vector12(daf::DAF, seg::SPKSegmentType20, time::Number)

    head = header(seg)
    data = cache(seg)

    # TODO: implement me

end


"""
    find_logical_record(head::SPKSegmentHeader2, time::Number)
"""
function find_logical_record(head::SPKSegmentHeader20, time::Number)
    
    # The index is returned in 0-based notation (i.e., the first record has index 0)
    index = floor(Int, (time - head.tstart)/head.tlen)

    if index == head.n 
        # This happens only when the time equals the final segment time
        index -= 1 
    end 

    return index
end

"""
    get_coefficients!(daf::DAF, head, cache, index::Int)
"""
function get_coefficients!(
    daf::DAF, head::SPKSegmentHeader20, cache::SPKSegmentCache20, index::Int
    )

    # Check whether the coefficients for this record are already loaded
    index == cache.id[1] && return nothing
    cache.id[1] = index 

    # Address of desired logical record 
    k = 8*(head.iaa-1) + head.recsize*index

    # For type 20 we do not have repeated midpoint and radius values 
    @inbounds for j = 1:3
        
        for i = 1:head.N 
            cache.A[j, i] = get_float(array(daf), k, endian(daf))
            k += 8
        end

        cache.p[j] = get_float(array(daf), k, endian(daf))
        k += 8

    end

    nothing 

end

function chebyshev_time(head::SPKSegmentHeader20, time::Number, index::Int)
    tbeg = head.tstart + head.tlen*index
    hlen = head.tlen/2 
    return (time - tbeg)/hlen - 1
end



function ∫cheb(cache::InterpCache, cₖ, t::Number, N::Int, Δl, Δt, tlen, p₀)

    # Retrieve the work buffer 
    t2 = 2t 

    if N >= 3 
        A₂ = cₖ[1, 1] - cₖ[1, 3]/2
    else 
        A₂ = cₖ[1, 1]
    end 

    adegp1 = N >= 3 ? cₖ[1, N-1]/(2N-2) : 0 
    adegp2 = N >= 2 ? cₖ[1, N]/(2N) : 0

    f3 = 0
    f2 = 0

    f1 = N == 1 ? A₂ : adegp2 

    z3 = 0 
    z2 = 0
    z1 = f1 

    w1 = 0 
    w2 = 0

    i = N
    while i > 1 

        if i == 2 
            Ai = A₂
        elseif i < N 
            Ai = (cₖ[1, i-1] - cₖ[1, i+1])/(2i-2)
        else 
            Ai = adegp1 
        end 

        f3 = f2 
        f2 = f1 
        f1 = Ai + t2*f2 - f3 

        z3 = z2 
        z2 = z1 
        z1 = Ai - z3 

        w3 = w2 
        w2 = w1 
        w1 = cₖ[1, i] + t2*w2 - w3

        i -= 1
    end

    c0 = z2 
    it = c0 + t*f1 - f2 
    
    println((t*f1 - f2)*tlen/2/Δt)

    vx = cₖ[1, 1] + t*w1 - w2
    x  = p₀[1] + it*tlen/2/Δt

    Δlt = Δl/Δt

    y, z = 0, 0 
    vy, vz = 0, 0

    return Δl*x, Δl*y, Δl*z, Δlt*vx, Δlt*vy, Δlt*vz

end 


function ∫chebyshev(cache::InterpCache, cₖ, t::Number, idx::Int, N::Int, Δl, Δt, tlen, p₀)

    # Retrieve the work buffer 
    Tₙ  = get_buffer(cache, 1, t)
    iTₙ = get_buffer(cache, 2, t)

    # idx is 0 or 3
    ix = 1 + idx 
    iy = 2 + idx 
    iz = 3 + idx 

    @inbounds begin 

        Tₙ[1] = 1
        Tₙ[2] = t 
        Tₙ[3] = 2t*t - 1

        iTₙ[1] = t 
        iTₙ[2] = t*t/2 

        vx = cₖ[ix, 1] + t*cₖ[ix, 2] + Tₙ[3]*cₖ[ix, 3]
        vy = cₖ[iy, 1] + t*cₖ[iy, 2] + Tₙ[3]*cₖ[iy, 3]
        vz = cₖ[iz, 1] + t*cₖ[iz, 2] + Tₙ[3]*cₖ[iz, 3]

        x = t*cₖ[ix, 1] + iTₙ[2]*cₖ[ix, 2] 
        y = t*cₖ[iy, 1] + iTₙ[2]*cₖ[iy, 2] 
        z = t*cₖ[iz, 1] + iTₙ[2]*cₖ[iz, 2]

        for j = 4:N

            Tₙ[j]    = 2*t*Tₙ[j-1] - Tₙ[j-2]
            iTₙ[j-1] = (Tₙ[j]/j - Tₙ[j-2]/(j-2))/2

            x += iTₙ[j-1]*cₖ[ix, j-1] 
            y += iTₙ[j-1]*cₖ[iy, j-1]
            z += iTₙ[j-1]*cₖ[iz, j-1]

            if j < N 
                vx += Tₙ[j]*cₖ[ix, j]
                vy += Tₙ[j]*cₖ[iy, j]
                vz += Tₙ[j]*cₖ[iz, j]
            end
        end

        println(x*tlen/Δt/2)

        # il /2 dovrebbe essere perché tu dividi per il radius ma qui tlen è 2*radius
        x = p₀[1] + x*tlen/Δt/2
        y = p₀[2] + y*tlen/Δt/2
        z = p₀[3] + z*tlen/Δt/2
        
    end

    Δlt = Δl/Δt

    return Δl*x, Δl*y, Δl*z, Δlt*vx, Δlt*vy, Δlt*vz

end 
