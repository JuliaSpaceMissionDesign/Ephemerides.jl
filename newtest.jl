include("src/Ephemeris.jl")

kernels = ["/home/michele/spice/kernels/spk/de440.bsp",
           "/home/michele/spice/kernels/spk/de430.bsp"]

# kernels = ["res/de421.bsp", "res/sat452.bsp"]

using BenchmarkTools
using FrameTransformations

using CALCEPH
using JSMDInterfaces.Ephemeris
using JSMDUtils.Math: D¹, D², D³
using CalcephEphemeris

eph = EphemerisProvider(kernels);
ephem = CalcephProvider(kernels[2]);

t1, t2, _ = ephem_spk_timespan(eph)

# TEST VECTOR 3
y1 = ephem_vector3(eph, 399, 3, 0);
y2 = @MVector zeros(3);
ephem_compute!(y2, ephem, DJ2000, 0, 3, 399, 0);
minimum(abs.(y1-y2))

# TEST VECTOR 6
y3 = ephem_vector6(eph, 399, 3, 0);
y4 = @MVector zeros(6);
ephem_compute!(y4, ephem, DJ2000, 0, 3, 399, 1);
minimum(abs.(y3-y4))

# TEST VECTOR 9
y5 = ephem_vector9(eph, 399, 3, 0);
y6 = @MVector zeros(9);
ephem_compute!(y6, ephem, DJ2000, 0, 3, 399, 2)
minimum(abs.(y5-y6))

# TEST VECTOR 12
y7 = ephem_vector12(eph, 399, 3, 0);
y8 = @MVector zeros(12);
ephem_compute!(y8, ephem, DJ2000, 0, 3, 399, 3)
minimum(abs.(y7-y8))

@benchmark ephem_vector3($eph, 399, 3, 0.0)
@benchmark ephem_compute!($y2, $ephem, $DJ2000, 0.0, 3, 399, 0)

@benchmark ephem_vector6($eph, 399, 3, 0.0)
@benchmark ephem_compute!($y4, $ephem, $DJ2000, 0.0, 3, 399, 1)

@benchmark ephem_vector9($eph, 399, 3, 0.0)
@benchmark ephem_compute!($y6, $ephem, $DJ2000, 0.0, 3, 399, 2)

@benchmark ephem_vector12($eph, 399, 3, 0.0)
@benchmark ephem_compute!($y8, $ephem, $DJ2000, 0.0, 3, 399, 3)








# D¹(t->spk_vector3(daf, seg, desc, t), 0) 
# D¹(t->spk_vector4(daf, seg, desc, t), 0) 

# D¹(t->spk_vector6(daf, seg, desc, t), 0) 

# spk_vector3(daf, seg, desc, 0.0)
# spk_vector6(daf, seg, desc, 0.0)

# spk_vector9(daf, seg, desc, 0.0)
# spk_vector12(daf, seg, desc, 0.0)

# function test_dual(x::Number)

#     SA[x, x^2, 3x^3 + 4]

# end

# function test_dual(x::Dual{T}) where T

#     pos = test_dual(value(x))
#     vel = SA[1, 2*value(x), 9*value(x)^2]

#     SA[
#         Dual{T}(pos[1], vel[1]*partials(x)), 
#         Dual{T}(pos[2], vel[2]*partials(x)), 
#         Dual{T}(pos[3], vel[3]*partials(x))
#     ]

# end

# D¹(t->test_dual(t), 1)
