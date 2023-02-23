
include("../spk/spk2.jl")

filename = "res/de421.bsp"

array = mmap(filename, Vector{UInt8});
header = load_daf(array)

summaries = get_daf_summaries(array, header);
for summ in summaries 
    println(parse_spk_segment_header(summ, header.lend))
end

spkhead = parse_spk_segment_header(summaries[1], header.lend)
segment = SPKSegmentType2(array, spkhead, header.lend)

# TESTING 
furnsh("/home/michele/spice/kernels/lsk/naif0012.tls", filename)


for i = 1:5000 
    et = rand(spkhead.tstart:spkhead.tend)

    p3 = spkpos(string(spkhead.tid), et, "J2000", "NONE", string(spkhead.cid))[1]

    @test p3 ≈ pos atol=1e-14 rtol=1e14

end

some_vector = [3.141592, 1.0, 2.0, 3.0]
B = pointer(some_vector)
recovered_vector = Base.unsafe_wrap(Array, B, 4)

ptr = Ptr{Float64}(pointer(array, 0))
unsafe_wrap(Vector{Float64}, ptr, 3)

reinterpret(reshape, Float64, array)

@views A = reinterpret(Float64, array[1:32])
A = unsafe_wrap(Array, Ptr{Float64}(ptr), (3, 10), own=false)


et = rand(spkhead.tstart:spkhead.tend)
@benchmark vector3($array, $header, $spkhead, $segment, $et)
@benchmark vector32($array, $header, $spkhead, $segment, $et)


