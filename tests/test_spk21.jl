
using Ephemerides
using Test 

using BenchmarkTools
using CalcephEphemeris
using Tempo
using StaticArrays

using JSMDInterfaces.Ephemeris
using JSMDUtils.Math: D¹, D², D³

kernel = "res/example1spk_seg21.bsp"

ephj = EphemerisProvider(kernel);
ephc = CalcephProvider(kernel);

# Center and target bodies 
cid = 10 
tid = 2025143

t1, t2 = ephem_spk_timespan(ephj)[1:2]
tj = (t1+t2)/2
tc = tj/86400

# Test values 
yc1 = @MVector zeros(3);
yc2 = @MVector zeros(6);

yj1 = ephem_vector3(ephj, cid, tid, tj);
yj2 = ephem_vector6(ephj, cid, tid, tj);

ephem_compute!(yc1, ephc, DJ2000, tc, tid, cid, 0);
ephem_compute!(yc2, ephc, DJ2000, tc, tid, cid, 1);

@test yj1 ≈ yc1 atol=1e-9 rtol=1e-9
@test yj2 ≈ yc2 atol=1e-9 rtol=1e-9

# Test performance 
@benchmark ephem_vector3($ephj, $cid, $tid, $tj)
@benchmark ephem_compute!($yc1, $ephc, $DJ2000, $tc, $tid, $cid, 0)

@benchmark ephem_vector6($ephj, $cid, $tid, $tj)
@benchmark ephem_compute!($yc2, $ephc, $DJ2000, $tc, $tid, $cid, 1)

# Test if AUTODIFF works 
D¹(t->ephem_vector3(ephj, cid, tid, t), tj)
D²(t->ephem_vector3(ephj, cid, tid, t), tj)
D³(t->ephem_vector3(ephj, cid, tid, t), tj)

D¹(t->ephem_vector6(ephj, cid, tid, t), tj)
D²(t->ephem_vector6(ephj, cid, tid, t), tj)
D³(t->ephem_vector6(ephj, cid, tid, t), tj)

# Test allocations with AD 
@benchmark D¹(t->ephem_vector3($ephj, $cid, $tid, t), $tj)
@benchmark D²(t->ephem_vector3($ephj, $cid, $tid, t), $tj)
@benchmark D³(t->ephem_vector3($ephj, $cid, $tid, t), $tj)

@benchmark D¹(t->ephem_vector6($ephj, $cid, $tid, t), $tj)
@benchmark D²(t->ephem_vector6($ephj, $cid, $tid, t), $tj)
@benchmark D³(t->ephem_vector6($ephj, $cid, $tid, t), $tj)


