"""
    AbstractEphemeris
Abstract type to represent ephemeris
"""
abstract type AbstractEphemeris end

"""
    SPKSegment
Abstract SPK Segment container.
"""
abstract type SPKSegment end

"""
    parse_segment(daf::DAF, name::String, symmary::Vector{UInt8}, ::Val)
SPK file segment parser.
"""
function parse_segment end 

function unsafe_position! end 
function unsafe_velocity! end 
function unsafe_state! end 

struct DAF
    filename::String
    array::Vector{UInt8}
    little::Bool
    id::String
    nd::Int32
    ni::Int32
    name::String
    first::Int32
    last::Int32
    ss::Int32
    nc::Int32
end

struct SPK <: AbstractEphemeris
    daf::DAF
    segments::Dict{Int,Dict{Int,SPKSegment}}
    graph::NodeGraph{Int, Int}
end

struct SPKDescriptor 
    name::String
    firstsec::Float64
    lastsec::Float64
    firstdate::Float64
    lastdate::Float64
    target::Int
    center::Int
    frame::Int
    spktype::Int
    firstaddr::Int
    lastaddr::Int
end

mutable struct Type02Segment <: SPKSegment
    descriptor::SPKDescriptor
    firstword::Int
    lastword::Int
    initialsecond::Float64
    intlen::Float64
    rsize::Int
    n_records::Int
    order::Int
    cached_record::Int
    cache::Matrix{Float64}
    x::Vector{Float64}
    t::Vector{Float64}
end
