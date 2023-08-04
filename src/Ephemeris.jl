
using Mmap

include("utils.jl")
include("daf.jl")

struct EphemerisProvider
    files::Vector{DAF}
    descriptors::Vector{DAFSegmentDescriptor}
end

EphemerisProvider(files::String) = EphemerisProvider([files])

function EphemerisProvider(files::Vector{String})

    # Initial parsing of each DAF file 
    nfiles = length(files)
    dafs = Vector{DAF}(undef, nfiles)

    seg_descr = DAFSegmentDescriptor[]

    @inbounds for j = nfiles:-1:1
        
        # Retrieve the header, comment, and content of each DAF 
        dafs[j] = DAF(files[j])

        # Retrieve all the summary records of this DAF  
        summaries = parse_daf_summaries(dafs[j])

        # Parse the descriptors 
        for k in eachindex(summaries)
            push!(seg_descr, DAFSegmentDescriptor(dafs[j], summaries[k], j))
        end

    end

    return EphemerisProvider(dafs, seg_descr)

end

kernels = ["/home/michele/spice/kernels/spk/de440.bsp",
           "/home/michele/spice/kernels/spk/de430.bsp"]

# using BenchmarkTools
eph = EphemerisProvider(kernels);

eph1 = EphemerisProvider(kernels[1]);
eph2 = EphemerisProvider(kernels[2]);

# @benchmark eph = EphemerisProvider($kernels)