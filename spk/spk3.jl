# TODO: missing test. Can you find a kernel of type 3?

# each segment can have different properties! 
struct SPKSegmentType3
    tstart::Float64         # segment starting epoch 
    tlen::Float64           # interval length 
    order::Int              # polynomial order (retrieved from array size)
    n::Int                  # number of records in segment 
    A::Matrix{Float64}
    x::Vector{Float64}
end

function SPKSegmentType3(array, spkhead::SPKSegmentHeader, lend::Bool)

    i0 = 8*(spkhead.faa-4)
    tstart = get_float(array, i0, lend)
    tlen = get_float(array, i0+8, lend)

    # polynomial order 
    rsize = Int(get_float(array, i0+16, lend))
    order = (rsize - 2) ÷ 6

    # number of records
    n = Int(get_float(array, i0+24, lend))

    # initialises arrays
    A = zeros(6, order)
    x = zeros(order)

    SPKSegmentType3(tstart, tlen, order, n, A, x)

end

function find_logical_record(time::Number, seg::SPKSegmentType3)
    # MISSING: Check that you are actually within the segment time boundaries! 
    idx, tfrac = divrem(time - seg.tstart, seg.tlen) 
    index = round(Int, idx)

    if index == seg.n 
        # This should only happen when time equals the final segment time
        index -= 1 
        tfrac = seg.tlen 
    end 

    return index, tfrac
end

function get_coefficients!(array, head::SPKSegmentHeader, seg::SPKSegmentType3, 
            index::Integer, lend::Bool)

    ncomp = 6 # components (3 for pos, 6 for vel)

    # size of each logical record in bytes for spk type 2
    recsize = 8*(seg.order*ncomp+2)

    # address of desired logical record (skipping mid and radius because they are all = )
    k = 8*(head.iaa-1) + recsize*index + 16

    @inbounds for j = 1:seg.order 
        for i = 1:ncomp 
            seg.A[i, j] = get_float(array, k + 8*(j-1) + 8*(i-1)*seg.order, lend)
        end
    end

end

function chebyshev!(seg::SPKSegmentType3, tfrac::Number)

    seg.x[1] = 1.0 
    seg.x[2] = 2*tfrac/seg.tlen - 1 
    @inbounds for i = 3:seg.order 
        seg.x[i] = 2.0*seg.x[2]*seg.x[i-1] - seg.x[i-2]
    end

end 

function vector3(array, daf::DAFHeader, head::SPKSegmentHeader, 
            seg::SPKSegmentType3, time::Number)

    index, tfrac = find_logical_record(time, seg)
    get_coefficients!(array, head, seg, index, daf.lend);
    chebyshev!(seg, tfrac)

    x, y, z = 0.0, 0.0, 0.0
    @inbounds @simd for i in eachindex(seg.x)
        x += seg.x[i]*seg.A[1, i]
        y += seg.x[i]*seg.A[2, i]
        z += seg.x[i]*seg.A[3, i]
    end

    return SA[x, y, z]

end

function vector6(array, daf::DAFHeader, head::SPKSegmentHeader, 
            seg::SPKSegmentType3, time::Number)
   
    index, tfrac = find_logical_record(time, seg)
    get_coefficients!(array, head, seg, index, daf.lend);
    chebyshev!(seg, tfrac)

    x, y, z = 0.0, 0.0, 0.0
    vx, vy, vz = 0.0, 0.0, 0.0

    @inbounds @simd for i in eachindex(seg.x)
        x += seg.x[i]*seg.A[1, i]
        y += seg.x[i]*seg.A[2, i]
        z += seg.x[i]*seg.A[3, i]
        vx += seg.x[i]*seg.A[4, i]
        vy += seg.x[i]*seg.A[5, i]
        vz += seg.x[i]*seg.A[6, i]
    end

    return SA[x, y, z, vx, vy, vz]
            
end