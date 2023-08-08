include("src/Ephemeris.jl")

kernels = ["/home/michele/spice/kernels/spk/de440.bsp",
           "/home/michele/spice/kernels/spk/de430.bsp"]


using BenchmarkTools
using FrameTransformations

using Tempo
using CALCEPH
using JSMDInterfaces.Ephemeris
using CalcephEphemeris

eph = EphemerisProvider(kernels);
ephem = CalcephProvider(kernels[2])


y1 = ephem_vector3(eph, 399, 3, 0.0)

y2 = @MVector zeros(3)
ephem_compute!(y2, ephem, DJ2000, 0.0, 3, 399, 0)


@benchmark ephem_vector3($eph, 399, 3, 0.0)
@benchmark ephem_compute!($y2, $ephem, $DJ2000, 0.0, 3, 399, 0)


# @code_warntype ephem_vector3(eph, 399, 3, 0.0)
# @code_warntype vector4(eph, 399, 3, 0.0)

# daf1 = eph.files[1];
# link = eph.spklinks[399][3][1]

# @code_warntype ephem_vector3(daf1, link, 0.0)
