
include("../spk/spk1.jl")

filename = "res/example1spk_seg1.bsp";   

array = mmap(filename, Vector{UInt8});
header = load_daf(array)

summaries = get_daf_summaries(array, header);
for summ in summaries 
    println(parse_spk_segment_header(summ, header.lend))
end

spkhead = parse_spk_segment_header(summaries[1], header.lend)
segment = SPKSegmentType1(array, spkhead, header.lend)

# TESTING 
furnsh("/home/michele/spice/kernels/lsk/naif0012.tls", filename)

for i = 1:5000 
    et = rand(spkhead.tstart:spkhead.tend)

    p3 = spkpos(string(spkhead.tid), et, "J2000", "NONE", string(spkhead.cid))[1]
    pos = vector3(array, header, spkhead, segment, et)

    @test p3 ≈ pos atol=1e-14 rtol=1e14

end

