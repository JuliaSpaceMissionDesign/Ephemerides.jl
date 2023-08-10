
include("src/Ephemeris.jl")

using BenchmarkTools
using CalcephEphemeris
using JSMDInterfaces.Ephemeris
using Tempo
using JSMDUtils.Math: D¹, D², D³
# using FrameTransformations
using Test

kernels = "res/example1spk_seg1.bsp";   

eph = EphemerisProvider(kernels);
ephem = CalcephProvider(kernels)

daf = get_daf(eph, 1)

seg = segment_list(daf).spk1[1]
desc = descriptors(daf)[1]

tms = initial_time(desc):final_time(desc)
for _ = 1:50000
    t = rand(tms)

    y1 = ephem_vector3(eph, 2000001, 0, t)
    y2 = @MVector zeros(3)
    ephem_compute!(y2, ephem, DJ2000, t/86400, 0, 2000001, 0)
    @test y1 ≈ y2 atol=1e-14 rtol=1e14
end

t = (6.669648e8 + 7.326288e8)/2
cid = Int(center(desc)) 
tid = Int(target(desc))

t = 0

y1 = ephem_vector3(eph, cid, tid, t)
y2 = @MVector zeros(3);
ephem_compute!(y2, ephem, DJ2000, t/86400, tid, cid, 0)
y1 - y2

@benchmark ephem_vector3($eph, $cid, $tid, $t)
@benchmark ephem_compute!($y2, $ephem, $DJ2000, $t/86400, $tid, $cid, 0)


D¹(t->ephem_vector3(eph, cid, tid, t), t)
D²(t->ephem_vector3(eph, cid, tid, t), t)

D¹(t->ephem_vector6(eph, cid, tid, t), t)
D²(t->ephem_vector6(eph, cid, tid, t), t)

@benchmark D¹(t->ephem_vector3($eph, $cid, $tid, t), $t)
@benchmark D²(t->ephem_vector3($eph, $cid, $tid, t), $t)

@benchmark D¹(t->ephem_vector6($eph, $cid, $tid, t), $t)
@benchmark D²(t->ephem_vector6($eph, $cid, $tid, t), $t)

y3 = ephem_vector6(eph, cid, tid, t)

y4 = @MVector zeros(6);
ephem_compute!(y4, ephem, DJ2000, t/86400, tid, cid, 1);
y4 - y3

@benchmark ephem_vector3($eph, $cid, $tid, $t)
@benchmark ephem_compute!($y4, $ephem, $DJ2000, $t/86400, $tid, $cid, 0)

@benchmark ephem_vector6($eph, $cid, $tid, $t)
@benchmark ephem_compute!($y4, $ephem, $DJ2000, $t/86400, $tid, $cid, 1)

@benchmark ephem_vector3($eph, 2000001, 0, 0.0)
@benchmark spk_vector3($daf, $seg, 0.0)

@code_warntype spk_vector3(daf, seg, 0.0)
@code_warntype find_logical_record(daf, header(seg), 0.0)
@code_warntype get_coefficients!(daf, header(seg), cache(seg), 0)
@code_warntype compute_mda3(cache(seg), 0)

@benchmark find_logical_record($daf, header($seg), 0.0)
@benchmark get_coefficients!($daf, header($seg), cache($seg), 0)
@benchmark compute_mda3(cache($seg), 0)



# spk_vector3(daf::DAF, seg::SPKSegmentType1, desc::DAFSegmentDescriptor, time::Number) 

# head = header(seg)
# cach = cache(seg)

# head.n 
# head.epochs

# # TESTING 
# furnsh("/home/michele/spice/kernels/lsk/naif0012.tls", filename)

# for i = 1:5000 
#     et = rand(spkhead.tstart:spkhead.tend)

#     p3 = spkpos(string(spkhead.tid), et, "J2000", "NONE", string(spkhead.cid))[1]
#     pos = vector3(array, header, spkhead, segment, et)

#     @test p3 ≈ pos atol=1e-14 rtol=1e14

# end

