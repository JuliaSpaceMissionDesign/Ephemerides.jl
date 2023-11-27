
# ----------------------------------
# ABSTRACT TYPES
# ----------------------------------

"""
    AbstractSPKSegment 

Abstract type for all SPK segment types.
"""
abstract type AbstractSPKSegment end 

"""
    spk_field(spk::AbstractSPKSegment)

Return the field number in the [`Ephemerides.SPKSegmentList`](@ref) associated to the given SPK 
segment type.
"""
function spk_field(::T) where {T <: AbstractSPKSegment}
    throw(ErrorException("`spk_field` must be implemented for SPK segment type $T"))
end

"""
    header(spk::AbstractSPKSegment)

Return the segment header.
"""
function header(::T) where {T <: AbstractSPKSegment}
    throw(ErrorException("`header` must be implemented for SPK segment type $T"))
end

"""
    cache(spk::AbstractSPKSegment)

Return the segment cache data.
"""
function cache(::T) where {T <: AbstractSPKSegment}
    throw(ErrorException("`cache` must be implemented for SPK segment type $T"))
end

"""
    AbstractSPKHeader 

Abstract type for all SPK segment type headers. 
"""
abstract type AbstractSPKHeader end 

""" 
    AbstractSPKCache 

Abstract type for all SPK segment type caches.
"""
abstract type AbstractSPKCache end 


# ----------------------------------
# SPK TYPE 1
# ----------------------------------

"""
    SPKSegmentHeader1 <: AbstractSPKHeader

Header instance for SPK segments of type 1 and 21.

### Fields 
- `n` -- `Int` number of records in the segment 
- `ndirs` -- `Int` number of directory epochs
- `epochs` -- Storage for directory epochs or epochs (when ndirs = 0)
- `iaa` - `Int` initial segment file address 
- `etid` -- `Int` initial address for the epoch table (after all the MDA records)
- `recsize` - `Int` Number of double numbers stored in each MDA record
- `maxdim` - `Int` MDA dimension (fixed to 15 for type 1)
"""
struct SPKSegmentHeader1 <: AbstractSPKHeader 
    n::Int                  
    ndirs::Int              
    epochs::Vector{Float64} 
    iaa::Int                
    etid::Int               
    recsize::Int            
    maxdim::Int            
end

""" 
    SPKSegmentCache1 <: AbstractSPKCache 

Cache instance for SPK segments of type 1 and 21. The fields contained within this cache 
are taken from the FORTRAN NAIF's SPICE implementation for type 1 SPK segments. 
"""
struct SPKSegmentCache1 <: AbstractSPKCache
    
    tl::Vector{Float64}
    g::Vector{Float64}
    refpos::Vector{Float64}
    refvel::Vector{Float64}
    dt::Matrix{Float64}
    kqmax::Vector{Int}
    kq::Vector{Int}
    id::MVector{1, Int}
    fc::DiffCache{Vector{Float64}, Vector{Float64}} 
    wc::DiffCache{Vector{Float64}, Vector{Float64}}     
    w::DiffCache{Vector{Float64},  Vector{Float64}}     
    vct::DiffCache{Vector{Float64}, Vector{Float64}}    

end

""" 
    SPKSegmentType1 <: AbstractSPKSegment

Segment instance for SPK segments of type 1 and 21, which contain Modified Difference Arrays 
(MDA). This data type is normally used for spacecraft whose ephemerides are produced by JPL's 
principal trajectory integrator DPTRAJ. 

### Fields 
- `head` -- Segment header 
- `cache` -- Segment cache 

### References 
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html)
"""
struct SPKSegmentType1 <: AbstractSPKSegment
    head::SPKSegmentHeader1
    cache::Vector{SPKSegmentCache1}
end

@inline header(spk::SPKSegmentType1) = spk.head 
@inline @inbounds cache(spk::SPKSegmentType1) = spk.cache[Threads.threadid()]


# ----------------------------------
# SPK TYPE 2
# ----------------------------------

"""
    SPKSegmentHeader2 <: AbstractSPKHeader

Header instance for SPK segments of type 2 and 3.

### Fields 
- `tstart` -- `Float64` initial epoch of the first record, in seconds since J2000
- `tlen` -- `Float64` interval length covered by each record, in seconds
- `order` -- `Int` polynomial order 
- `N` -- `Int` number of coefficients in each window
- `n` -- `Int` number of records in the segment
- `recsize` -- `Int` byte size of each logical record
- `ncomp` -- `Int` number of vector components
- `iaa` -- `Int` initial segment file address
- `type` -- `Int` SPK segment type, either 2 or 2
"""
struct SPKSegmentHeader2 <: AbstractSPKHeader
    tstart::Float64     
    tlen::Float64       
    order::Int
    N::Int
    n::Int              
    recsize::Int
    ncomp::Int
    iaa::Int    
    type::Int
end

""" 
    SPKSegmentCache2 <: AbstractSPKCache 

Cache instance for SPK segments of type 2 and 3.

### Fields 
- `A` -- Chebyshev's polynomial coefficients, with size (ncomp, order)
- `p` -- Stores the record mid point and radius and scale factor
- `buff` -- Stores the buffers for the Chebyshev polynomials
- `id` -- Index of the currently loaded logical record
"""
struct SPKSegmentCache2 <: AbstractSPKCache
    A::Matrix{Float64}
    p::MVector{3, Float64} 
    buff::InterpCache{Float64}
    id::MVector{1, Int}
end

""" 
    SPKSegmentType2 <: AbstractSPKSegment

Segment instance for SPK segments of type 2 and 3, which contain Chebyshev polynomial 
coefficients for the position and/or state of the body as function of time. This data type 
is normally used for planet barycenters, and for satellites whose ephemerides are integrated.

### Fields 
- `head` -- Segment header 
- `cache` -- Segment cache 

### References 
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html)
"""
struct SPKSegmentType2 <: AbstractSPKSegment
    head::SPKSegmentHeader2
    cache::Vector{SPKSegmentCache2}
end

@inline header(spk::SPKSegmentType2) = spk.head 
@inline @inbounds cache(spk::SPKSegmentType2) = spk.cache[Threads.threadid()]


# ----------------------------------
# SPK TYPE 5
# ----------------------------------

"""
    SPKSegmentHeader5 <: AbstractSPKHeader

Header instance for SPK segments of type 5.

### Fields 

"""
struct SPKSegmentHeader5 <: AbstractSPKHeader
    GM::Float64 
    n::Int 
    ndirs::Int 
    etid::Int
    epochs::Vector{Float64}
end

""" 
    SPKSegmentCache5 <: AbstractSPKCache 

Cache instance for SPK segments of type 5.

### Fields 
- `id` -- Index of the currently loaded logical record
"""
struct SPKSegmentCache5 <: AbstractSPKCache
    s1::MVector{6, Float64}
    s2::MVector{6, Float64}
    ts::MVector{2, Float64}
    id::MVector{1, Int}
end

""" 
    SPKSegmentType5 <: AbstractSPKSegment

Segment instance for SPK segments of type 5. 

### Fields 
- `head` -- Segment header 
- `cache` -- Segment cache 

### References 
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html)
"""
struct SPKSegmentType5 <: AbstractSPKSegment
    head::SPKSegmentHeader5
    cache::Vector{SPKSegmentCache5}
end

@inline header(spk::SPKSegmentType5) = spk.head 
@inline @inbounds cache(spk::SPKSegmentType5) = spk.cache[Threads.threadid()]


# ----------------------------------
# SPK TYPE 8
# ----------------------------------

"""
    SPKSegmentHeader8 <: AbstractSPKHeader

Header instance for SPK segments of type 8 and 12.

### Fields 
- `tstart` -- `Float64` segment starting epoch, in TDB seconds since J2000 
- `tlen` -- `Float64` interval length, in seconds
- `order` -- `Int` interpolating polynomial degree
- `N` -- `Int` group size
- `n` -- `Int` number of states in the segment
- `iaa` - `Int` initial segment file address 
- `iseven` -- `Bool` true for even group size
- `type` -- `Int` SPK type (either 8 or 12)
"""
struct SPKSegmentHeader8 <: AbstractSPKHeader
    tstart::Float64     
    tlen::Float64       
    order::Int
    N::Int
    n::Int
    iaa::Int      
    iseven::Bool
    type::Int
end

"""
    SPKSegmentCache8 <: AbstractSPKCache

Cache instance for SPK segments of type 8 and 12.
"""
struct SPKSegmentCache8 <: AbstractSPKCache
    states::Matrix{Float64}
    buff::InterpCache{Float64}
    id::MVector{1, Int}
end 

""" 
    SPKSegmentType8 <: AbstractSPKSegment

Segment instance for SPK segments of type 8 and 12.

### Fields 
- `head` -- Segment header 
- `cache` -- Segment cache 

### References 
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html)
"""
struct SPKSegmentType8 <: AbstractSPKSegment
    head::SPKSegmentHeader8
    cache::Vector{SPKSegmentCache8}
end

@inline header(spk::SPKSegmentType8) = spk.head 
@inline @inbounds cache(spk::SPKSegmentType8) = spk.cache[Threads.threadid()]


# ----------------------------------
# SPK TYPE 9
# ----------------------------------

"""
    SPKSegmentHeader9 <: AbstractSPKHeader

Header instance for SPK segments of type 9 and 13.

### Fields 
- `n` -- `Int` number of states in the segment
- `ndirs` -- `Int` number of epoch directories
- `epochs` -- Storage for directory epochs or epochs (when ndirs = 0)
- `iaa` - `Int` initial segment file address 
- `etid` -- `Int` initial address for the epoch table (after all the state data)
- `order` -- `Int` interpolating polynomial degree
- `N` -- `Int` group size 
- `iseven` -- `Bool` true for even group size
- `type` -- `Int` SPK type (either 9 or 13)
"""
struct SPKSegmentHeader9 <: AbstractSPKHeader
    n::Int
    ndirs::Int 
    epochs::Vector{Float64}
    iaa::Int
    etid::Int 
    order::Int 
    N::Int 
    iseven::Bool 
    type::Int
end

"""
    SPKSegmentCache9 <: AbstractSPKCache

Cache instance for SPK segments of type 9 and 13.
"""
struct SPKSegmentCache9 <: AbstractSPKCache
    epochs::Vector{Float64}
    states::Matrix{Float64}
    buff::InterpCache{Float64}
    id::MVector{1, Int}
end 

""" 
    SPKSegmentType9 <: AbstractSPKSegment

Segment instance for SPK segments of type 9 and 13.

### Fields 
- `head` -- Segment header 
- `cache` -- Segment cache 

### References 
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html)
"""
struct SPKSegmentType9 <: AbstractSPKSegment
    head::SPKSegmentHeader9
    cache::Vector{SPKSegmentCache9}
end

@inline header(spk::SPKSegmentType9) = spk.head 
@inline @inbounds cache(spk::SPKSegmentType9) = spk.cache[Threads.threadid()]


# ----------------------------------
# SPK TYPE 14
# ----------------------------------

"""
    SPKSegmentHeader14 <: AbstractSPKHeader

Header instance for SPK segments of type 14.

### Fields 
- `order` -- `Int` interpolating polynomial degree
- `n` -- `Int` number of packets in the segment
- `ndirs` -- `Int` number of epoch directories
- `epochs` -- Storage for directory epochs or epochs (when ndirs = 0)
- `etid` -- `Int` initial address for the epoch table (after all the state data)
- `ptid` -- `Int` initial address for the packet table (after the constants)
- `pktsize` -- `Int` size of each data packet excluding the packet information area.
- `pktoff` -- `Int` offset of the packet data from the packet start 
- `ncomp` -- `Int` number of states coefficients (= 6 for SPK 14)
- `N` -- `Int` number of polynomial coefficients
"""
struct SPKSegmentHeader14 <: AbstractSPKHeader
    order::Int
    n::Int 
    ndirs::Int 
    epochs::Vector{Float64}
    etid::Int 
    ptid::Int  
    pktsize::Int 
    pktoff::Int  
    ncomp::Int
    N::Int 
end

""" 
    SPKSegmentType14 <: AbstractSPKSegment

Segment instance for SPK segments of type 14.

### Fields 
- `head` -- Segment header 
- `cache` -- Segment cache 

### References 
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html)
"""
struct SPKSegmentType14 <: AbstractSPKSegment 
    head::SPKSegmentHeader14 
    cache::Vector{SPKSegmentCache2}
end

@inline header(spk::SPKSegmentType14) = spk.head 
@inline @inbounds cache(spk::SPKSegmentType14) = spk.cache[Threads.threadid()]


# ----------------------------------
# SPK TYPE 17
# ----------------------------------

"""
    SPKSegmentHeader17 <: AbstractSPKHeader

Header instance for SPK segments of type 17.

### Fields 
"""
struct SPKSegmentHeader17 <: AbstractSPKHeader
    epoch::Float64 
    sma::Float64 
    h::Float64 
    k::Float64 
    lon::Float64
    p::Float64 
    q::Float64 
    dlpdt::Float64
    dmldt::Float64 
    dnodedt::Float64
    ra::Float64
    de::Float64
    R::SMatrix{3, 3, Float64, 9}
end


""" 
    SPKSegmentType17 <: AbstractSPKSegment

Segment instance for SPK segments of type 17.

### Fields 
- `head` -- Segment header 

!!! note 
    SPK segments of type 17 do not require a cache because they do not extract any 
    additional coefficients at runtime.

### References 
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html)
"""
struct SPKSegmentType17 <: AbstractSPKSegment
    head::SPKSegmentHeader17
end

@inline header(spk::SPKSegmentType17) = spk.head 


# ----------------------------------
# SPK TYPE 18
# ----------------------------------

"""
    SPKSegmentHeader18 <: AbstractSPKHeader

Header instance for SPK segments of type 18.

### Fields 
- `n` -- `Int` number of states in the segment
- `ndirs` -- `Int` number of epoch directories
- `epochs` -- Storage for directory epochs or epochs (when ndirs = 0)
- `iaa` - `Int` initial segment file address 
- `etid` -- `Int` initial address for the epoch table (after all the state data)
- `order` -- `Int` interpolating polynomial degree
- `N` -- `Int` group size 
- `subtype` -- `Int` type 18 subtype, either 0 (Hermite) or 1 (Lagrange)
- `packetsize` -- `Int` packet size for each point, either 12 (Hermite) or 6 (Lagrange)
"""
mutable struct SPKSegmentHeader18 <: AbstractSPKHeader
    n::Int
    ndirs::Int 
    epochs::Vector{Float64}
    iaa::Int
    etid::Int 
    order::Int 
    N::Int 
    subtype::Int
    packetsize::Int
end

"""
    SPKSegmentCache18 <: AbstractSPKCache

Cache instance for SPK segments of type 18.
"""
struct SPKSegmentCache18 <: AbstractSPKCache
    p::MVector{3, Int}
    epochs::Vector{Float64}
    states::Matrix{Float64}
    buff::InterpCache{Float64}
end 


# ----------------------------------
# SPK TYPE 19
# ----------------------------------

"""
    SPKSegmentHeader19 <: AbstractSPKHeader

Header instance for SPK segments of type 19.

### Fields 
- `n` -- `Int` number of states in the segment.
- `ndirs` -- `Int` number of epoch directories.
- `times` -- Storage for interval directories or start times (when ndirs = 0).
- `iaa` - `Int` initial segment file address.
- `etid` -- `Int` byte address for the interval table (after all the minisegment data).
- `ptid` -- `Int` byte for the pointer table.
- `usefirst` -- `Bool` boundary flag, true if the preceding segment should be used.
- `type` -- `Int` either type 18 or 19.
"""
struct SPKSegmentHeader19 <: AbstractSPKHeader
    n::Int 
    ndirs::Int 
    times::Vector{Float64}
    iaa::Int 
    etid::Int 
    ptid::Int 
    usefirst::Bool
    type::Int 
end

"""
    SPKSegmentCache19 <: AbstractSPKCache

Cache instance for SPK segments of type 19.
"""
struct SPKSegmentCache19 <: AbstractSPKCache
    minihead::SPKSegmentHeader18
    minidata::SPKSegmentCache18
    id::MVector{1, Int}
end 

""" 
    SPKSegmentType19 <: AbstractSPKSegment

Segment instance for SPK segments of type 18 and 19. Type 18 segments are treated as 
special cases of a type 19 with a single mini-segment.

### Fields 
- `head` -- Segment header 
- `cache` -- Segment cache 

### References 
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html)
"""
struct SPKSegmentType19 <: AbstractSPKSegment
    head::SPKSegmentHeader19
    cache::Vector{SPKSegmentCache19}
end

@inline header(spk::SPKSegmentType19) = spk.head 
@inline @inbounds cache(spk::SPKSegmentType19) = spk.cache[Threads.threadid()]


# ----------------------------------
# SPK TYPE 20
# ----------------------------------

"""
    SPKSegmentHeader20 <: AbstractSPKHeader

Header instance for SPK segments of type 20.

### Fields 
- `dscale` -- `Float64` length conversion factor
- `tscale` -- `Float64` time conversion factor
- `tstart` -- `Float64` initial epoch of the first record, in seconds since J2000
- `tlen` -- `Float64` interval length covered by each record, in seconds
- `recsize` -- `Int` byte size of each logical record
- `order` -- `Int` polynomial order 
- `N` -- `Int` number of coefficients in each window
- `n` -- `Int` number of records in the segment
- `iaa` -- `Int` initial segment file address
"""
struct SPKSegmentHeader20 <: AbstractSPKHeader 
    dscale::Float64 
    tscale::Float64 
    tstart::Float64
    tlen::Float64 
    recsize::Int 
    order::Int 
    N::Int 
    n::Int 
    iaa::Int 
end

""" 
    SPKSegmentCache2 <: AbstractSPKCache 

Cache instance for SPK segments of type 20.

### Fields 
- `id` -- Index of the currently loaded logical record
- `p` -- Stores the record position constants
- `A` -- Chebyshev's polynomial coefficients, with size (ncomp, order)
- `buff` -- Stores the buffers for the Chebyshev polynomials

"""
struct SPKSegmentCache20 <: AbstractSPKCache
    id::MVector{1, Int}
    p::MVector{3, Float64}
    A::Matrix{Float64}
    buff::InterpCache{Float64}
end

""" 
    SPKSegmentType2 <: AbstractSPKSegment

Segment instance for SPK segments of type 20, which contain Chebyshev polynomial coefficients 
for the position and/or state of the body as function of time. This data type is normally 
used for planet barycenters, and for satellites whose ephemerides are integrated.

### Fields 
- `head` -- Segment header 
- `cache` -- Segment cache 

### References 
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html)
"""
struct SPKSegmentType20 <: AbstractSPKSegment
    head::SPKSegmentHeader20
    cache::Vector{SPKSegmentCache20}
end

@inline header(spk::SPKSegmentType20) = spk.head 
@inline @inbounds cache(spk::SPKSegmentType20) = spk.cache[Threads.threadid()]


"""
    SPK_SEGMENT_MAPPING

A dictionary mapping SPK segment types to the field index of the [`SPKSegmentList`](@ref).
"""
const SPK_SEGMENTLIST_MAPPING = Dict(
    1 => 1,
    2 => 2,
    3 => 2,
    8 => 3,
    9 => 4,
    12 => 3,
    13 => 4,
    14 => 5,
    18 => 6,
    19 => 6,
    20 => 7,
    21 => 1,
    17 => 8,
    5 => 9
)

# ----------------------------------
# SPK SEGMENT LIST
# ----------------------------------

"""
    SPKSegmentList 

A container object to efficiently store all the different SPK segments that are contained 
within a single DAF file.

---

    SPKSegmentList()

Initialises an empty `SPKSegmentList` object.

### See also 
See also [`Ephemerides.add_segment!`](@ref)

"""
struct SPKSegmentList
    
    spk1::Vector{SPKSegmentType1}
    spk2::Vector{SPKSegmentType2}
    spk8::Vector{SPKSegmentType8}
    spk9::Vector{SPKSegmentType9}
    spk14::Vector{SPKSegmentType14}
    spk19::Vector{SPKSegmentType19}
    spk20::Vector{SPKSegmentType20}
    spk17::Vector{SPKSegmentType17}
    spk5::Vector{SPKSegmentType5}

    function SPKSegmentList()
        new(
            SPKSegmentType1[], 
            SPKSegmentType2[], 
            SPKSegmentType8[],
            SPKSegmentType9[],
            SPKSegmentType14[],
            SPKSegmentType19[],
            SPKSegmentType20[],
            SPKSegmentType17[], 
            SPKSegmentType5[],
        )
    end
end

"""
    add_segment!(list::SPKSegmentList, spk::AbstractSPKSegment)

Add the SPK segment to the proper vector within the given [`Ephemerides.SPKSegmentList`](@ref) `list` 
"""
@inline function add_segment!(list::SPKSegmentList, spk::AbstractSPKSegment)
    push!(getfield(list, spk_field(spk)), spk)
    nothing
end

"""
    get_segment(list::SPKSegmentList, lid::Int, eid::Int)

Return the segment contained in the `lid` list at index `eid`.
"""
@inline get_segment(list::SPKSegmentList, lid::Int, eid::Int) = getfield(list, lid)[eid]
