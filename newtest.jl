include("src/Ephemeris.jl")

kernels = ["/home/michele/spice/kernels/spk/de440.bsp",
           "/home/michele/spice/kernels/spk/de430.bsp"]


using BenchmarkTools
using FrameTransformations

using CALCEPH
using JSMDInterfaces.Ephemeris
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


# @code_warntype ephem_vector3(eph, 399, 3, 0.0)
# @code_warntype vector4(eph, 399, 3, 0.0)

# daf1 = eph.files[1];
# link = eph.spklinks[399][3][1]

# @code_warntype ephem_vector3(daf1, link, 0.0)
