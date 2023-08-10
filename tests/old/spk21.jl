
using CalcephEphemeris
using SPICE
using Test
using Tempo 
using JSMDInterfaces.Ephemeris

include("../spk/reader.jl")
include("../spk/spk21.jl")

filename = "res/example1spk_seg21.bsp"

array = mmap(filename, Vector{UInt8});
header = load_daf(array)

summaries = get_daf_summaries(array, header);
for summ in summaries 
    println(parse_spk_segment_header(summ, header.lend))
end

spkhead = parse_spk_segment_header(summaries[1], header.lend)
segment = SPKSegmentType21(array, spkhead, header.lend)

# TESTING 
furnsh("/home/michele/spice/kernels/lsk/naif0012.tls", filename)

ephem = CalcephProvider(filename)
p4 = @MVector zeros(3)

for i = 1:5000 
    et = rand(spkhead.tstart:spkhead.tend)

    p3 = spkpos(string(spkhead.tid), et, "J2000", "NONE", string(spkhead.cid))[1]
    pos = vector3(array, header, spkhead, segment, et)
    
    ephem_compute!(p4, ephem, DJ2000, et/86400, Int(spkhead.tid), Int(spkhead.cid), 0)

    @test p3 ≈ pos atol=1e-14 rtol=1e14
    # println(maximum(abs.(p3 - pos)))
    println(maximum(abs.(p4 - pos)))

end

