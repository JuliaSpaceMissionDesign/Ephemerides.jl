
# ----------------------------------
# ABSTRACT TYPES
# ----------------------------------

"""
    AbstractSPKSegment 

Abstract type for all SPK segment types.
"""
abstract type AbstractSPKSegment end 

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

struct SPKSegmentHeader1 <: AbstractSPKHeader 
    n::Int      # number of records in the segment
    ndirs::Int 
end

struct SPKSegmentCache1 <: AbstractSPKCache
    
    tl::Vector{Float64}
    g::Vector{Float64}
    refpos::Vector{Float64}
    refvel::Vector{Float64}
    dt::Matrix{Float64}
    kqmax::Vector{Int}
    kq::Vector{Int}
    fc::Vector{Float64}
    wc::Vector{Float64}
    w::Vector{Float64}

end

struct SPKSegmentType1
    head::SPKSegmentHeader1
    cache::SPKSegmentCache1
end


# ----------------------------------
# SPK TYPE 2
# ----------------------------------
# TODO: add number of components to be able to handle TYPE 3 aswell with same header! 
struct SPKSegmentHeader2 <: AbstractSPKHeader
    tstart::Float64     # segment starting epoch 
    tlen::Float64       # interval length 
    order::Int          # polynomial order (retrieved from array size)
    n::Int              # number of records in segment 
end

struct SPKSegmentCache2 <: AbstractSPKCache
    A::Matrix{Float64}
    x::Vector{Float64}
end

struct SPKSegmentType2 <: AbstractSPKSegment
    head::SPKSegmentHeader2
    cache::SPKSegmentCache2
end



# ----------------------------------
# SPK SEGMENT LIST
# ----------------------------------

struct SPKSegmentList
    
    spk1::Vector{SPKSegmentType1}
    spk2::Vector{SPKSegmentType2}

    function SPKSegmentList()
        new(
            SPKSegmentType1[], 
            SPKSegmentType2[]
        )
    end
end
