

# each segment can have different properties! 
struct SPKSegmentType8
    tstart::Float64         # segment starting epoch 
    tlen::Float64           # interval length 
    order::Int              # polynomial degree
    N::Int                  # group size (order + 1)
    n::Int                  # number of states in segment 
    states::Matrix{Float64}
    work::Vector{Float64}
end

function SPKSegmentType8(array, spkhead::SPKSegmentHeader, lend::Bool)

    i0 = 8*(spkhead.faa-4)

    tstart = get_float(array, i0, lend)
    tlen = get_float(array, i0+8, lend)

    # polynomial order 
    order = Int(get_float(array, i0 + 16, lend))
    N = order + 1 # number of states

    # Number of states within the record
    n = Int(get_float(array, i0 + 24, lend))

    # initialises arrays
    states = zeros(N, 6)
    work = zeros(N)

    SPKSegmentType8(tstart, tlen, order, N, n, states, work)

end

function find_logical_record(time::Number, seg::SPKSegmentType8)
    # MISSING: Check that you are actually within the segment time boundaries! 
    
    Δt = time - seg.tstart 

    if seg.N % 2 == 0 # even group size 
        low = Int(Δt ÷ seg.tlen) + 1
        first = low - (seg.N ÷ 2) + 1
    else # odd group size  
        near = round(Int, Δt/seg.tlen) + 1
        first = near - seg.order ÷ 2
    end

    first_idx = min(max(1, first), seg.n-seg.order)

    return first_idx
end

function get_coefficients!(array, head::SPKSegmentHeader, seg::SPKSegmentType8, 
            index::Integer, lend::Bool)

    
    i0 = 8*(head.iaa - 1) + 8*6*(index-1)

    @inbounds for j = 1:6
        for i = 1:seg.N 
            seg.states[i, j] = get_float(array, i0 + 48*(i-1) + 8*(j-1), lend)
        end
    end

end

function lagrange(seg::SPKSegmentType8, newx::Number, istate::Integer)

    # This function is valid only for equally-spaced polynomials!

    # newx is a re-work of the ascissa that starts at 1
    # istate is the index of the desired state
    
    @inbounds for i = 1:seg.N 
        seg.work[i] = seg.states[i, istate]
    end

    # compute the lagrange polynomials using a recursive relationship
    @inbounds for j = 1:seg.N-1 
        for i = 1:seg.N-j 
            c1 = i + j - newx
            c2 = newx - i 

            seg.work[i] = (c1*seg.work[i] + c2*seg.work[i+1])/j
        end
    end

    return seg.work[1]

end 

function vector3(array, daf::DAFHeader, head::SPKSegmentHeader, 
            seg::SPKSegmentType8, time::Number)

    index = find_logical_record(time, seg)
    get_coefficients!(array, head, seg, index, daf.lend);

    tbegin = seg.tstart + (index-1)*seg.tlen 
    newx = (time - tbegin)/seg.tlen + 1

    x = lagrange(seg, newx, 1)
    y = lagrange(seg, newx, 2)
    z = lagrange(seg, newx, 3)

    return SA[x, y, z]

end

function vector6(array, daf::DAFHeader, head::SPKSegmentHeader, 
            seg::SPKSegmentType8, time::Number)

    index = find_logical_record(time, seg)
    get_coefficients!(array, head, seg, index, daf.lend);

    tbegin = seg.tstart + (index-1)*seg.tlen 
    newx = (time - tbegin)/seg.tlen + 1

    x = lagrange(seg, newx, 1)
    y = lagrange(seg, newx, 2)
    z = lagrange(seg, newx, 3)

    vx = lagrange(seg, newx, 4)
    vy = lagrange(seg, newx, 5)
    vz = lagrange(seg, newx, 6)

    return SA[x, y, z, vx, vy, vz]

end