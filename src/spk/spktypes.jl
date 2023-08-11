
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
- `n` -- `Int` number of records in the segment
- `recsize` -- `Int` byte size of each logical record
- `ncomp` -- `Int` number of vector components
- `scale` -- `Float64` scale factor for the Chebyshev derivatives
- `iaa` -- `Int` initial segment file address
"""
struct SPKSegmentHeader2 <: AbstractSPKHeader
    tstart::Float64     
    tlen::Float64       
    order::Int
    n::Int              
    recsize::Int
    ncomp::Int
    scale::Float64
    iaa::Int    
end

""" 
    SPKSegmentCache2 <: AbstractSPKCache 

Cache instance for SPK segments of type 2 and 3.

### Fields 
- `A` -- Chebyshev's polynomial coefficients, with size (ncomp, order)
- `x1` -- Values of the Chebyshev's polynomials
- `x2` -- Derivatives of the Chebyshev's polynomials
- `id` -- Index of the currently loaded logical record
"""
struct SPKSegmentCache2 <: AbstractSPKCache
    A::Matrix{Float64} 
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

Header instance for SPK segments of type 8.

### Fields 
- `tstart` -- `Float64` segment starting epoch, in TDB seconds since J2000 
- `tlen` -- `Float64` interval length, in seconds
- `order` -- `Int` interpolating polynomial degree
- `N` -- `Int` group size (order + 1)
- `n` -- `Int` number of states in the segment
- `iaa` - `Int` initial segment file address 
- `iseven` -- `Bool` true for even group size
"""
struct SPKSegmentHeader8 <: AbstractSPKHeader
    tstart::Float64     
    tlen::Float64       
    order::Int
    N::Int
    n::Int
    iaa::Int      
    iseven::Bool
end

"""
    SPKSegmentCache8 <: AbstractSPKCache

Cache instance for SPK segments of type 8.
"""
struct SPKSegmentCache8 <: AbstractSPKCache
    states::Matrix{Float64}
    work::Vector{Float64}
    id::MVector{1, Int}
end 

""" 
    SPKSegmentType8 <: AbstractSPKSegment

Segment instance for SPK segments of type 8.

### Fields 
- `head` -- Segment header 
- `cache` -- Segment cache 

### References 
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html)
"""
struct SPKSegmentType8 <: AbstractSPKSegment
    header::SPKSegmentHeader8
    cache::Vector{SPKSegmentCache8}
end

@inline header(spk::SPKSegmentType8) = spk.header 
@inline @inbounds cache(spk::SPKSegmentType8) = spk.cache[Threads.threadid()]


"""
    SPK_SEGMENT_MAPPING

A dictionary mapping SPK segment types to the field index of the [`SPKSegmentList`](@ref).
"""
const SPK_SEGMENTLIST_MAPPING = Dict(
    1 => 1,
    2 => 2,
    3 => 3,
    21 => 1,
    8 => 4
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

    function SPKSegmentList()
        new(
            SPKSegmentType1[], 
            SPKSegmentType2[], 
            SPKSegmentType3[], 
            SPKSegmentType8[]
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
