include("src/Ephemeris.jl")

kernels = ["/home/michele/spice/kernels/spk/de440.bsp",
           "/home/michele/spice/kernels/spk/de430.bsp"]

# using BenchmarkTools
eph = EphemerisProvider(kernels);

# eph1 = EphemerisProvider(kernels[1]);
# eph2 = EphemerisProvider(kernels[2]);

using BenchmarkTools
using FrameTransformations

using CALCEPH
using JSMDInterfaces.Ephemeris
using CalcephEphemeris

ephem = CalcephProvider(kernels[2])
# prefetch(ephem)

# compute(eph::Ephem,jd0::Float64,time::Float64,
   # target::Integer,center::Integer,unit::Integer)

ephem_vector3(eph, 399, 3, 0.0)
# compute(ephem, Float64(DJ2000), 0.0, 3, 399, useNaifId+unitKM+unitSec, 0)

y = @MVector zeros(3)
ephem_compute!(y, ephem, DJ2000, 0.0, 3, 399, 0)

# jd0 = Float64(DJ2000)
# unit = useNaifId+unitKM+unitSec

# using BenchmarkTools
@benchmark ephem_vector3($eph, 399, 3, 0.0)
@benchmark ephem_compute!($y, $ephem, $DJ2000, 0.0, 3, 399, 0)


@code_warntype ephem_vector3(eph, 399, 3, 0.0)
@code_warntype vector4(eph, 399, 3, 0.0)

daf1 = eph.files[1];
link = eph.spklinks[399][3][1]

@code_warntype ephem_vector3(daf1, link, 0.0)
