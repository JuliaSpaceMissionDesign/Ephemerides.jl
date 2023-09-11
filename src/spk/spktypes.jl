
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

Return the field number in the [`SPKSegmentList`](@ref) associated to the given SPK 
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
"""
struct SPKSegmentHeader1 <: AbstractSPKHeader 
    n::Int                  # number of records in the segment
    ndirs::Int              # Number of directory epochs
    epochs::Vector{Float64} # Vector storing directory epochs or epochs themselves (when ndirs = 0)
    iaa::Int                # Initial segment address
    etid::Int               # Initial address for the epoch table (after all the records)
    recsize::Int            # Number of double numbers stored in each MDA record
    maxdim::Int             # Dimensions 
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
end

""" 
    SPKSegmentCache2 <: AbstractSPKCache 

Cache instance for SPK segments of type 2 and 3.

### Fields 
- `A` -- Chebyshev's polynomial coefficients, with size (ncomp, order)
- `p` -- Stores the record mid point and radius and scale factor
- `x1` -- Values of the Chebyshev's polynomials
- `x2` -- Derivatives of the Chebyshev's polynomials
- `id` -- Index of the currently loaded logical record
"""
struct SPKSegmentCache2 <: AbstractSPKCache
    A::Matrix{Float64}
    p::MVector{3, Float64} 
    x1::DiffCache{Vector{Float64}, Vector{Float64}}
    x2::DiffCache{Vector{Float64}, Vector{Float64}}
    id::MVector{1, Int}
end

""" 
    SPKSegmentType2 <: AbstractSPKSegment

Segment instance for SPK segments of type 2, which contain Chebyshev polynomial coefficients 
for the position and/or state of the body as function of time. This data type is normally 
used for planet barycenters, and for satellites whose ephemerides are integrated.

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
# SPK TYPE 3
# ----------------------------------

"""
    SPKSegmentType3 <: AbstractSPKSegment
"""
struct SPKSegmentType3 <: AbstractSPKSegment
    head::SPKSegmentHeader2 
    cache::Vector{SPKSegmentCache2}
end

@inline header(spk::SPKSegmentType3) = spk.head 
@inline @inbounds cache(spk::SPKSegmentType3) = spk.cache[Threads.threadid()]


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
- `N` -- `Int` group size (order + 1)
- `n` -- `Int` number of states in the segment
- `iaa` - `Int` initial segment file address 
- `iseven` -- `Bool` true for even group size
- `type` -- `Int` SPK type
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
# SPK TYPE 18
# ----------------------------------

"""
    SPKSegmentHeader18 <: AbstractSPKHeader

Header instance for SPK segments of type 18.

### Fields 

"""
struct SPKSegmentHeader18 <: AbstractSPKHeader
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

""" 
    SPKSegmentType18 <: AbstractSPKSegment

Segment instance for SPK segments of type 18.

### Fields 
- `head` -- Segment header 
- `cache` -- Segment cache 

### References 
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html)
"""
struct SPKSegmentType18 <: AbstractSPKSegment
    head::SPKSegmentHeader18
    cache::Vector{SPKSegmentCache18}
end

@inline header(spk::SPKSegmentType18) = spk.head 
@inline @inbounds cache(spk::SPKSegmentType18) = spk.cache[Threads.threadid()]


"""
    SPK_SEGMENT_MAPPING

A dictionary mapping SPK segment types to the field index of the [`SPKSegmentList`](@ref).
"""
const SPK_SEGMENTLIST_MAPPING = Dict(
    1 => 1,
    2 => 2,
    3 => 3,
    8 => 4,
    9 => 5,
    12 => 4,
    13 => 5,
    18 => 6,
    21 => 1,
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
See also [`add_segment!`](@ref)

"""
struct SPKSegmentList
    
    spk1::Vector{SPKSegmentType1}
    spk2::Vector{SPKSegmentType2}
    spk3::Vector{SPKSegmentType3}
    spk8::Vector{SPKSegmentType8}
    spk9::Vector{SPKSegmentType9}
    spk18::Vector{SPKSegmentType18}

    function SPKSegmentList()
        new(
            SPKSegmentType1[], 
            SPKSegmentType2[], 
            SPKSegmentType3[], 
            SPKSegmentType8[],
            SPKSegmentType9[],
            SPKSegmentType18[]
        )
    end
end

"""
    add_segment!(list::SPKSegmentList, spk::AbstractSPKSegment)

Add the SPK segment to the proper vector within the given [`SPKSegmentList`](@ref) `list` 
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
