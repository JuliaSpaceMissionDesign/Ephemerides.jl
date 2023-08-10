

include("spk/spk2.jl")

# each segment can have different properties! 
struct SPKSegmentDualType2
    tstart::Float64 # segment starting epoch 
    tlen::Float64   # interval length 
    order::Int      # polynomial order (retrieved from array size)
    n::Int          # number of records in segment 
    A::Matrix{Float64}
    x::DiffCache{Vector{Float64}, Vector{Float64}}
end

function SPKSegmentDualType2(array, spkhead::SPKSegmentHeader, lend::Bool)

    i0 = 8*(spkhead.faa-4)
    tstart = get_float(array, i0, lend)
    tlen = get_float(array, i0+8, lend)

    # polynomial order 
    rsize = Int(get_float(array, i0+16, lend))

    #  the order of the polynomial is actually = (order-1)
    order = (rsize - 2) ÷ 3

    # number of records
    n = Int(get_float(array, i0+24, lend))

    # initialises arrays
    A = zeros(3, order)
    x = zeros(order)

    SPKSegmentDualType2(tstart, tlen, order, n, A, DiffCache(x))

end

function find_logical_record(time::Number, seg::SPKSegmentDualType2)
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

function get_coefficients!(array, head::SPKSegmentHeader, seg::SPKSegmentDualType2, 
            index::Integer, lend::Bool)

    ncomp = 3 # components (3 for pos, 6 for vel)

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

function chebyshev(seg::SPKSegmentDualType2, time::Number, tfrac::Number)

    x = get_tmp(seg.x, time)

    x[1] = 1.0 
    x[2] = 2*tfrac/seg.tlen - 1 
    @inbounds for i = 3:seg.order 
        x[i] = 2.0*x[2]*x[i-1] - x[i-2]
    end

    return x

end 

@views function vector133(array, daf::DAFHeader, head::SPKSegmentHeader, 
            seg::SPKSegmentDualType2, time::T) where {T}

    index, tfrac = find_logical_record(time, seg)
    get_coefficients!(array, head, seg, index, daf.lend);

    xseg = chebyshev(seg, time, tfrac)

    x, y, z = T(0.0), T(0.0), T(0.0)
    @inbounds @simd for i in eachindex(xseg)
        x += xseg[i]*seg.A[1, i]
        y += xseg[i]*seg.A[2, i]
        z += xseg[i]*seg.A[3, i]
    end

    return SA[x, y, z]
end


filename = "res/de421.bsp"

array = mmap(filename, Vector{UInt8});
header = load_daf(array)

summaries = get_daf_summaries(array, header);
for summ in summaries 
    println(parse_spk_segment_header(summ, header.lend))
end

spkhead = parse_spk_segment_header(summaries[1], header.lend);
segment = SPKSegmentDualType2(array, spkhead, header.lend);

# TESTING 
furnsh("/home/michele/spice/kernels/lsk/naif0012.tls", filename)

et = rand(spkhead.tstart:spkhead.tend)

pos = vector133(array, header, spkhead, segment, et)

vector133(array, header, spkhead, segment, et)
D¹(t->vector133(array, header, spkhead, segment, t), et)
D²(t->vector133(array, header, spkhead, segment, t), et)
D³(t->vector133(array, header, spkhead, segment, t), et)

@benchmark vector133($array, $header, $spkhead, $segment, $et)
@benchmark D¹(t->vector133($array, $header, $spkhead, $segment, t), $et)
@benchmark D²(t->vector133($array, $header, $spkhead, $segment, t), $et)
@benchmark D³(t->vector133($array, $header, $spkhead, $segment, t), $et)



