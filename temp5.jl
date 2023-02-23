using Test

using BenchmarkTools
using SPICE

# NOTE: fatalerror("order>=2 is not supported on segment of type 5");

include("spk/reader.jl")

filename = "res/example1spk_seg5.bsp"

array = mmap(filename, Vector{UInt8});
header = load_daf(array)

summaries = get_daf_summaries(array, header);
for summ in summaries 
    println(parse_spk_segment_header(summ, header.lend))
end

spkhead = parse_spk_segment_header(summaries[1], header.lend)
# segment = SPKSegmentType5(array, spkhead, header.lend)

et = 0.5*(spkhead.tend + spkhead.tstart)

# epoch x y z vx vy vz

# find epochs that bracket the desired time et 

eph = Ephem(filename)
pr = positionRecords(eph)

pos = compute(eph, pr[1].startEpoch, 0.0, pr[1].target, pr[1].center, unitKM+unitSec+useNaifId, 1)

# TESTING 
furnsh("/home/michele/spice/kernels/lsk/naif0012.tls", filename)


for i = 1:5000 
    et = rand(spkhead.tstart:spkhead.tend)

    p3 = spkezr(string(spkhead.tid), et, "J2000", "NONE", string(spkhead.cid))[1]
    pos = vector6(array, header, spkhead, segment, et)

    @test p3 ≈ pos atol=1e-14 rtol=1e14

end
