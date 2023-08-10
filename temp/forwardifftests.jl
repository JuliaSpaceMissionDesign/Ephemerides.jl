include("src/Ephemeris.jl")

using BenchmarkTools
using FrameTransformations

using CALCEPH
using JSMDInterfaces.Ephemeris
using JSMDUtils.Math: D¹, D², D³
using CalcephEphemeris

kernels = ["/home/michele/spice/kernels/spk/de440.bsp"]
eph = EphemerisProvider(kernels);

et = 0.0
@benchmark D¹(t->ephem_vector3($eph, 399, 3, t), $et)
@benchmark D²(t->ephem_vector3($eph, 399, 3, t), $et)
@benchmark D³(t->ephem_vector3($eph, 399, 3, t), $et)

@benchmark D¹(t->ephem_vector6($eph, 399, 3, t), $et)
@benchmark D²(t->ephem_vector6($eph, 399, 3, t), $et)
@benchmark D³(t->ephem_vector6($eph, 399, 3, t), $et)

@benchmark D¹(t->ephem_vector9($eph, 399, 3, t), $et)
@benchmark D²(t->ephem_vector9($eph, 399, 3, t), $et)
@benchmark D³(t->ephem_vector9($eph, 399, 3, t), $et)

@benchmark D¹(t->ephem_vector12($eph, 399, 3, t), $et)
@benchmark D²(t->ephem_vector12($eph, 399, 3, t), $et)
@benchmark D³(t->ephem_vector12($eph, 399, 3, t), $et)


@code_warntype ephem_vector3(eph, 399, 3, 0.0)
@code_warntype ephem_vector6(eph, 399, 3, 0.0)
@code_warntype ephem_vector9(eph, 399, 3, 0.0)

@code_warntype D¹(t->ephem_vector9(eph, 399, 3, t), 0.0)

D¹(t->ephem_vector3(eph, 399, 3, t), 0.0) - y7[4:6]
D²(t->ephem_vector3(eph, 399, 3, t), 0.0) - y7[7:9]
D³(t->ephem_vector3(eph, 399, 3, t), 0.0) - y7[10:12]

D¹(t->ephem_vector6(eph, 399, 3, t), 0.0) - y7[4:9]
D²(t->ephem_vector6(eph, 399, 3, t), 0.0) - y7[7:12]
D³(t->ephem_vector6(eph, 399, 3, t), 0.0) 

D¹(t->ephem_vector9(eph, 399, 3, t), 0.0) 
D²(t->ephem_vector9(eph, 399, 3, t), 0.0) 
D³(t->ephem_vector9(eph, 399, 3, t), 0.0) 

D¹(t->ephem_vector12(eph, 399, 3, t), 0.0) 
D²(t->ephem_vector12(eph, 399, 3, t), 0.0) 
D³(t->ephem_vector12(eph, 399, 3, t), 0.0) 


from = 399
to = 3
link = spk_links(eph)[399][3][1]

t = 0.0

daf = get_daf(eph, file_id(link));
seg = get_segment(segment_list(daf), 2, element_id(link));
desc = descriptor(link);

@code_warntype spk_vector3(daf, seg, desc, t)
@code_warntype spk_vector6(daf, seg, desc, t)
@code_warntype spk_vector9(daf, seg, desc, t)
@code_warntype spk_vector12(daf, seg, desc, t)

D¹(t->spk_vector3(daf, seg, desc, t), 0.0)
D²(t->spk_vector3(daf, seg, desc, t), 0.0)
D³(t->spk_vector3(daf, seg, desc, t), 0.0)

D¹(t->spk_vector6(daf, seg, desc, t), 0.0) 
D²(t->spk_vector6(daf, seg, desc, t), 0.0) 
D³(t->spk_vector6(daf, seg, desc, t), 0.0) 

D¹(t->spk_vector9(daf, seg, desc, t), 0.0) 
D²(t->spk_vector9(daf, seg, desc, t), 0.0) 
D³(t->spk_vector9(daf, seg, desc, t), 0.0) 

D¹(t->spk_vector12(daf, seg, desc, t), 0.0) 
D²(t->spk_vector12(daf, seg, desc, t), 0.0) 
D³(t->spk_vector12(daf, seg, desc, t), 0.0) 


@code_warntype D¹(t->spk_vector3(daf, seg, desc, t), 0.0)
@code_warntype D¹(t->spk_vector6(daf, seg, desc, t), 0.0) 
@code_warntype D¹(t->spk_vector9(daf, seg, desc, t), 0.0) 

@code_warntype D¹(t->spk_vector6(daf, seg, desc, t), 0.0) 
@code_warntype D²(t->spk_vector6(daf, seg, desc, t), 0.0) 
@code_warntype D³(t->spk_vector6(daf, seg, desc, t), 0.0) 

@code_warntype D¹(t->spk_vector9(daf, seg, desc, t), 0.0) 
@code_warntype D²(t->spk_vector9(daf, seg, desc, t), 0.0) 
@code_warntype D³(t->spk_vector9(daf, seg, desc, t), 0.0) 

@code_warntype D¹(t->spk_vector12(daf, seg, desc, t), 0.0) 
@code_warntype D²(t->spk_vector12(daf, seg, desc, t), 0.0) 
@code_warntype D³(t->spk_vector12(daf, seg, desc, t), 0.0) 