
include("src/Ephemeris.jl")
include("src/spk/spk3.jl")

using BenchmarkTools
using FrameTransformations

using CALCEPH
using JSMDInterfaces.Ephemeris
using JSMDUtils.Math: D¹, D², D³
using CalcephEphemeris

kernels = ["/home/michele/spice/kernels/spk/de440.bsp"]
eph = EphemerisProvider(kernels);

daf = get_daf(eph, 1)
desc = descriptors(daf)[10]

desc = DAFSegmentDescriptor(3, desc.tstart, desc.tend, desc.tid, desc.cid, 
        desc.axesid, desc.iaa, desc.faa)

seg = SPKSegmentType3(daf, desc)

t = 0.0
spk_vector3(daf, seg, t)
spk_vector6(daf, seg, t)
spk_vector9(daf, seg, t)
spk_vector12(daf, seg, t)

@code_warntype spk_vector3(daf, seg, t)
@code_warntype spk_vector6(daf, seg, t)
@code_warntype spk_vector9(daf, seg, t)
@code_warntype spk_vector12(daf, seg, t)

# spk_vector9(daf, seg, t)
# spk_vector12(daf, seg, t)

@benchmark D¹(t->spk_vector3($daf, $seg, t), 0.0)
@benchmark D²(t->spk_vector3($daf, $seg, t), 0.0)
@benchmark D³(t->spk_vector3($daf, $seg, t), 0.0)

@benchmark D¹(t->spk_vector6(daf, seg, t), 0.0)
@benchmark D²(t->spk_vector6(daf, seg, t), 0.0)
@benchmark D³(t->spk_vector6(daf, seg, t), 0.0)

@benchmark D¹(t->spk_vector9($daf, $seg, t), 0.0)
@benchmark D²(t->spk_vector9($daf, $seg, t), 0.0)
@benchmark D³(t->spk_vector9($daf, $seg, t), 0.0)


# D¹(t->spk_vector9(daf, seg, desc, t), 0.0)
# D¹(t->spk_vector12(daf, seg, desc, t), 0.0)
