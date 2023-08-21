
using Ephemerides
using Test 

using BenchmarkTools
using CalcephEphemeris
using Tempo
using StaticArrays

using JSMDInterfaces.Ephemeris
using JSMDUtils.Math: D¹, D², D³

kernel = "res/example1spk_seg8.bsp"

ephj = EphemerisProvider(kernel);
ephc = CalcephProvider(kernel);

ts, te = ephem_timespan(ephj)[1:2]
tj = (ts+te)/2

# Test values 
yc1 = @MVector zeros(3);
yc2 = @MVector zeros(6);
yc3 = @MVector zeros(9);
yc4 = @MVector zeros(12);

yj1 = ephem_vector3(ephj, 5, 0, tj);
yj2 = ephem_vector6(ephj, 5, 0, tj);
yj3 = ephem_vector9(ephj, 5, 0, tj);
yj4 = ephem_vector12(ephj, 5, 0, tj);

ephem_compute!(yc1, ephc, DJ2000, tj/86400, 0, 5, 0);
ephem_compute!(yc2, ephc, DJ2000, tj/86400, 0, 5, 1);
ephem_compute!(yc3, ephc, DJ2000, tj/86400, 0, 5, 2);
ephem_compute!(yc4, ephc, DJ2000, tj/86400, 0, 5, 3);
 
@test yj1 ≈ yc1 atol=1e-9 rtol=1e-9
@test yj2 ≈ yc2 atol=1e-9 rtol=1e-9
@test yj3 ≈ yc3 atol=1e-9 rtol=1e-9
@test yj4 ≈ yc4 atol=1e-9 rtol=1e-9

# Test performance 
@benchmark ephem_vector3($ephj, 5, 0, $tj)
@benchmark ephem_compute!($yc1, $ephc, $DJ2000, $tj/86400, 0, 5, 0)

@benchmark ephem_vector6($ephj, 5, 0, $tj)
@benchmark ephem_compute!($yc2, $ephc, $DJ2000, $tj/86400, 0, 5, 1)

@benchmark ephem_vector9($ephj, 5, 0, $tj)
@benchmark ephem_compute!($yc3, $ephc, $DJ2000, $tj/86400, 0, 5, 2)

@benchmark ephem_vector12($ephj, 5, 0, $tj)
@benchmark ephem_compute!($yc4, $ephc, $DJ2000, $tj/86400, 0, 5, 3)

# Test if AUTODIFF works 
D¹(t->ephem_vector3(ephj, 5, 0, t), tj)
D²(t->ephem_vector3(ephj, 5, 0, t), tj)
D³(t->ephem_vector3(ephj, 5, 0, t), tj)

D¹(t->ephem_vector6(ephj, 5, 0, t), tj)
D²(t->ephem_vector6(ephj, 5, 0, t), tj)
D³(t->ephem_vector6(ephj, 5, 0, t), tj)

D¹(t->ephem_vector9(ephj, 5, 0, t), tj)
D²(t->ephem_vector9(ephj, 5, 0, t), tj)
D³(t->ephem_vector9(ephj, 5, 0, t), tj)

D¹(t->ephem_vector12(ephj, 5, 0, t), tj)
D²(t->ephem_vector12(ephj, 5, 0, t), tj)
D³(t->ephem_vector12(ephj, 5, 0, t), tj)

# Test allocations with AD 
@benchmark D¹(t->ephem_vector3($ephj, 5, 0, t), $tj)
@benchmark D²(t->ephem_vector3($ephj, 5, 0, t), $tj)
@benchmark D³(t->ephem_vector3($ephj, 5, 0, t), $tj)

@benchmark D¹(t->ephem_vector6($ephj, 5, 0, t), $tj)
@benchmark D²(t->ephem_vector6($ephj, 5, 0, t), $tj)
@benchmark D³(t->ephem_vector6($ephj, 5, 0, t), $tj)

@benchmark D¹(t->ephem_vector9($ephj, 5, 0, t), $tj)
@benchmark D²(t->ephem_vector9($ephj, 5, 0, t), $tj)
@benchmark D³(t->ephem_vector9($ephj, 5, 0, t), $tj)

@benchmark D¹(t->ephem_vector12($ephj, 5, 0, t), $tj)
@benchmark D²(t->ephem_vector12($ephj, 5, 0, t), $tj)
@benchmark D³(t->ephem_vector12($ephj, 5, 0, t), $tj)


daf = Ephemerides.get_daf(ephj, 1)
link = ephj.spklinks[0][5][1]
seg = getfield(daf.seglist, link.lid)[link.eid]

x1 = Ephemerides.spk_vector4(daf, seg, 0.0)
x2 = Ephemerides.spk_vector3(daf, seg, 0.0)

@benchmark Ephemerides.spk_vector4($daf, $seg, 0.0)
@benchmark Ephemerides.spk_vector3($daf, $seg, 0.0)