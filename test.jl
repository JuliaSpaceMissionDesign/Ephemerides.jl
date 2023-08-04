

macro ephemeris(files)

    println(files)
    return quote 
        $files
    end
end


kernels = @ephemeris "/home/michele/spice/kernels/spk/de440.bsp"

kernels = @ephemeris begin 
    "/home/michele/spice/kernels/spk/de440.bsp",
    "/home/michele/spice/kernels/pck/moon_pa_de440_200625.bpc"
end

using CALCEPH
using FrameTransformations

eph1 = Ephem("/home/michele/spice/kernels/spk/de440.bsp")
eph2 = Ephem("/home/michele/spice/kernels/spk/de430.bsp")
eph3 = Ephem([
    "/home/michele/spice/kernels/spk/de430.bsp", 
    "/home/michele/spice/kernels/spk/de440.bsp"]
)

eph4 = Ephem([
    "/home/michele/spice/kernels/spk/de440.bsp", 
    "/home/michele/spice/kernels/spk/de430.bsp"]
)

a = compute(eph1, Float64(DJ2000), 0.0, 399, 301, useNaifId+unitKM+unitSec, 0)
b = compute(eph2, Float64(DJ2000), 0.0, 399, 301, useNaifId+unitKM+unitSec, 0)

c = compute(eph3, Float64(DJ2000), 0.0, 399, 301, useNaifId+unitKM+unitSec, 0)
d = compute(eph4, Float64(DJ2000), 0.0, 399, 301, useNaifId+unitKM+unitSec, 0)

using SPICE

furnsh("/home/michele/spice/kernels/spk/de440.bsp")
e = spkpos("399", DJ2000, "J2000", "NONE", "301")[1]
kclear()

furnsh("/home/michele/spice/kernels/spk/de430.bsp")
f = spkpos("399", DJ2000, "J2000", "NONE", "301")[1]
kclear()

furnsh("/home/michele/spice/kernels/spk/de440.bsp", "/home/michele/spice/kernels/spk/de430.bsp")
g = spkpos("399", DJ2000, "J2000", "NONE", "301")[1]
kclear()

furnsh("/home/michele/spice/kernels/spk/de430.bsp", "/home/michele/spice/kernels/spk/de440.bsp")
h = spkpos("399", DJ2000, "J2000", "NONE", "301")[1]
kclear()


macro ephemeris(files)
    
    if isa(files, String)
        println("single file")
    else
        # println(length(files.args[2].args))
        println(files.args[2].args[1][1] == "/")
    end
    
    # println(files.args[2].args)
    # println(length(files.args[2].args))
    # println(files[1])

    file = eval(files)
    # println(files.head)
    
    return quote 
         $files 
    end

end
 

eph = @ephemeris begin 
    "/home/michele/spice/kernels/spk/de440.bsp",
    "/home/michele/spice/kernels/pck/moon_pa_de440_200625.bpc"
end;

eph = @ephemeris "/home/michele/spice/kernels/spk/de440.bsp"


Meta.@dump begin 
    "/home/michele/spice/kernels/spk/de440.bsp",
    "/home/michele/spice/kernels/pck/moon_pa_de440_200625.bpc"
end

Meta.@dump "/home/michele/spice/kernels/spk/de440.bsp"
