
export EphemerisProvider

"""
    EphemerisProvider(file::String)
    EphemerisProvider(files::Vector{String})

Create an `EphemerisProvider` instance by loading a single or multiple binary ephemeris 
kernel files specified by `files`. Currently, only NAIF Double precision Array File (DAF)
kernels (i.e., SPK and PCK) are accepted.

### Example 
```julia-repl 
julia> eph = EphemerisProvider("PATH_TO_KERNEL")
EphemerisProvider([...])

julia> eph = EphemerisProvider(["PATH_TO_KERNEL_1", "PATH_TO_KERNEL_2"])
EphemerisProvider([])
```
"""
struct EphemerisProvider <: jEph.AbstractEphemerisProvider
    files::Vector{DAF}
    spklinks::SPKLinkTable
    pcklinks::SPKLinkTable
end

EphemerisProvider(files::AbstractString) = EphemerisProvider([files])
# TODO: warning if no file is loaded/ the kernel is empty?
function EphemerisProvider(files::Vector{<:AbstractString})

    # Initial parsing of each DAF file 
    ndafs = length(files)
    
    dafs = Vector{DAF}(undef, ndafs)
    @inbounds for fid = ndafs:-1:1
        dafs[fid] = DAF(files[fid])
        initialise_segments!(dafs[fid])
    end

    spklinks, pcklinks = create_linktables(dafs)
    return EphemerisProvider(dafs, spklinks, pcklinks) 

end

"""
    daf(eph::EphemerisProvider)

Return the [`DAF`](@ref) files stored in the ephemeris provider. 
"""
@inline get_daf(eph::EphemerisProvider) = eph.files

"""
    daf(eph::EphemerisProvider, id::Int)

Return the [`DAF`](@ref) file in the ephemeris provider at index `id`.
"""
@inline get_daf(eph::EphemerisProvider, id::Int) = get_daf(eph)[id]

"""
    spk_links(eph::EphemerisProvider)

Return the [`SPKLinkTable`] for the SPK segments.
"""
@inline spk_links(eph::EphemerisProvider) = eph.spklinks

"""
    pck_links(eph::EphemerisProvider)

Return the [`SPKLinkTable`](@ref) for the PCK segments.
"""
@inline pck_links(eph::EphemerisProvider) = eph.pcklinks