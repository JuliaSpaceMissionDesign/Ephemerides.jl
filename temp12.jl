using Test

using BenchmarkTools
using SPICE

include("spk/reader.jl")
filename = "/res/example1spk_seg12.bsp"

array = mmap(filename, Vector{UInt8});
header = load_daf(array)

summaries = get_daf_summaries(array, header);
for summ in summaries 
    println(parse_spk_segment_header(summ, header.lend))
end

spkhead = parse_spk_segment_header(summaries[1], header.lend)
# segment = SPKSegmentType8(array, spkhead, header.lend)

et = 0.5*(spkhead.tend + spkhead.tstart)
i0 = 8*(spkhead.faa-4)

tstart = get_float(array, i0, header.lend)
stepsize = get_float(array, i0 + 8, header.lend)
wndsize = Int(get_float(array, i0 + 16, header.lend)) + 1
degree = 2*wndsize - 1

N = degree + 1 # group size!

nstates = Int(get_float(array, i0 + 24, header.lend))

# Find first and start epochs of the states to be extracted!
Δt = et - tstart

if N % 2 == 0 # group size pari ! 
    low = Int(Δt ÷ stepsize) + 1
    first = low - (N ÷ 2) + 1
else # dispari !
    near = round(Int, Δt/stepsize) + 1
    first = near - degree/2
end

first_idx = min(max(1, first), nstates-degree)
last_idx = first_idx + degree

# Extract the states 
states = zeros(N, 6);

i0 = 8*(spkhead.iaa - 1) + 8*6*(first_idx-1)
for j = 1:6
    for i = 1:N
        states[i, j] = get_float(array, i0 + 8*6*(i-1) + 8*(j-1), header.lend)
    end
end

# compute polynomials
work = zeros(2N, 2)
posb = zeros(6)

# newx = (et - tstart) ÷ stepsize + 1
# newx = 1
# newx = 1
tbegin = tstart + (first_idx-1)*stepsize
newx = (et - tbegin)/stepsize + 1
for k = 1:3

    for i = 1:2:2N-1 
        work[i, 1] = states[i, k]
    end

    for i = 2:2:2N 
        work[i, 1] = states[i, k]*tlen
    end

    for i = 1:N-1 

        c1 = (i+1) - newx 
        c2 = newx - i 

        prev = 2i - 1
        this = prev + 1 
        next = this + 1

        work[prev, 2] = work[this, 1]
        work[this, 2] = work[next, 1] - work[prev, 1]

        temp = c2*work[this, 1] + work[prev, 1]
        
        work[this, 1] = c1*work[prev, 1] + c2*work[next, 1]
        work[prev, 1] = temp

    end

    work[2N-1, 2] = work[2N, 1]
    work[2N-1, 1] = work[2N, 1]*(newx-N) + work[2N-1, 1]

    for j = 2:2N-1 
        for i = 1:2N-j 
            xi = (i + 1)/2 
            xij = (i + j + 1) / 2 

            c1 = xij - newx 
            c2 = newx - xi 

            denom = xij - xi 

            work[i, 2] = (c1*work[i, 2] + c2*work[i+1, 2] + (work[i+1, 1] - work[i, 1]))/denom
            work[i, 1] = (c1*work[i, 1] + c2*work[i+1, 1]) / denom 
        end
    end

    posb[k] = work[1, 1]
    posb[k+3] = work[1, 2] / tlen

end

# TESTING 
furnsh("/home/michele/spice/kernels/lsk/naif0012.tls", filename)

et = spkhead.tstart
p3 = spkezr(string(spkhead.tid), et, "J2000", "NONE", string(spkhead.cid))[1]
pos = vector6(array, header, spkhead, segment, et)
p3-pos


for i = 1:5000 
    et = rand(spkhead.tstart:spkhead.tend)

    p3 = spkezr(string(spkhead.tid), et, "J2000", "NONE", string(spkhead.cid))[1]
    pos = vector6(array, header, spkhead, segment, et)

    @test p3 ≈ pos atol=1e-14 rtol=1e14

end
