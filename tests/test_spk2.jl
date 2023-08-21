
using Ephemerides
using Test 

using BenchmarkTools
using CalcephEphemeris
using Tempo
using StaticArrays

using JSMDInterfaces.Ephemeris
using JSMDUtils.Math: D¹, D², D³

kernel = "res/de421.bsp"

ephj = EphemerisProvider(kernel);
ephc = CalcephProvider(kernel);

# Test values 
yc1 = @MVector zeros(3);
yc2 = @MVector zeros(6);
yc3 = @MVector zeros(9);
yc4 = @MVector zeros(12);

yj1 = ephem_vector3(ephj, 399, 3, 0);
yj2 = ephem_vector6(ephj, 399, 3, 0);
yj3 = ephem_vector9(ephj, 399, 3, 0);
yj4 = ephem_vector12(ephj, 399, 3, 0);

ephem_compute!(yc1, ephc, DJ2000, 0, 3, 399, 0);
ephem_compute!(yc2, ephc, DJ2000, 0, 3, 399, 1);
ephem_compute!(yc3, ephc, DJ2000, 0, 3, 399, 2);
ephem_compute!(yc4, ephc, DJ2000, 0, 3, 399, 3);
 
@test yj1 ≈ yc1 atol=1e-9 rtol=1e-9
@test yj2 ≈ yc2 atol=1e-9 rtol=1e-9
@test yj3 ≈ yc3 atol=1e-9 rtol=1e-9
@test yj4 ≈ yc4 atol=1e-9 rtol=1e-9

# Test performance 
@benchmark ephem_vector3($ephj, 399, 3, 0.0)
@benchmark ephem_compute!($yc1, $ephc, $DJ2000, 0.0, 3, 399, 0)

@benchmark ephem_vector6($ephj, 399, 3, 0.0)
@benchmark ephem_compute!($yc2, $ephc, $DJ2000, 0.0, 3, 399, 1)

@benchmark ephem_vector9($ephj, 399, 3, 0.0)
@benchmark ephem_compute!($yc3, $ephc, $DJ2000, 0.0, 3, 399, 2)

@benchmark ephem_vector12($ephj, 399, 3, 0.0)
@benchmark ephem_compute!($yc4, $ephc, $DJ2000, 0.0, 3, 399, 3)

# Test if AUTODIFF works 
D¹(t->ephem_vector3(ephj, 399, 3, t), 0.0)
D²(t->ephem_vector3(ephj, 399, 3, t), 0.0)
D³(t->ephem_vector3(ephj, 399, 3, t), 0.0)

D¹(t->ephem_vector6(ephj, 399, 3, t), 0.0)
D²(t->ephem_vector6(ephj, 399, 3, t), 0.0)
D³(t->ephem_vector6(ephj, 399, 3, t), 0.0)

D¹(t->ephem_vector9(ephj, 399, 3, t), 0.0)
D²(t->ephem_vector9(ephj, 399, 3, t), 0.0)
D³(t->ephem_vector9(ephj, 399, 3, t), 0.0)

D¹(t->ephem_vector12(ephj, 399, 3, t), 0.0)
D²(t->ephem_vector12(ephj, 399, 3, t), 0.0)
D³(t->ephem_vector12(ephj, 399, 3, t), 0.0)

# Test allocations with AD 
@benchmark D¹(t->ephem_vector3($ephj, 399, 3, t), 0.0)
@benchmark D²(t->ephem_vector3($ephj, 399, 3, t), 0.0)
@benchmark D³(t->ephem_vector3($ephj, 399, 3, t), 0.0)

@benchmark D¹(t->ephem_vector6($ephj, 399, 3, t), 0.0)
@benchmark D²(t->ephem_vector6($ephj, 399, 3, t), 0.0)
@benchmark D³(t->ephem_vector6($ephj, 399, 3, t), 0.0)

@benchmark D¹(t->ephem_vector9($ephj, 399, 3, t), 0.0)
@benchmark D²(t->ephem_vector9($ephj, 399, 3, t), 0.0)
@benchmark D³(t->ephem_vector9($ephj, 399, 3, t), 0.0)

@benchmark D¹(t->ephem_vector12($ephj, 399, 3, t), 0.0)
@benchmark D²(t->ephem_vector12($ephj, 399, 3, t), 0.0)
@benchmark D³(t->ephem_vector12($ephj, 399, 3, t), 0.0)

