using LinearAlgebra: transpose!

const SIZE_FLOAT64 = sizeof(Float64)
const SIZE_INT32 = sizeof(Int32)

jd(sec) = 2451545.0 + sec / 86400.0
seconds(jd) = (jd - 2451545.0) * 86400.0

struct OutOfRangeError <: Exception
    date::Float64
    startdate::Float64
    finaldate::Float64
end

Base.showerror(io::IO, err::OutOfRangeError) = print(io,
   "The requested date $(err.date) is outside the interval ($(err.startdate), $(err.finaldate)).")

function get_segments_naifids(segments::Dict{Int, Dict{Int, SPKSegment}})
    s = Vector{Vector{Int}}()
    for (k,v) in segments
        for l in keys(v)
            push!(s, [k, l])
        end
    end
    sort!(s)
end

@inline function checkdate(seg::SPKSegment, tdb::Float64)
    if !(seg.descriptor.firstdate <= tdb <= seg.descriptor.lastdate)
        throw(OutOfRangeError(tdb, seg.descriptor.firstdate, seg.descriptor.lastdate))
    end
end

function SPK(filename)
    daf = DAF(filename)
    segments = Dict{Int, Dict{Int, SPKSegment}}()
    for (name, summary) in get_summaries(daf)
        spktype = reinterpret_getindex(Int32, summary, 29, daf.little)
        # TODO : important to add support for type 21 
        # TODO : extend support to type 3 and one of the simple types (12, 13, or)
        seg = parse_segment(daf, name, summary, Val(2))
        if haskey(segments, seg.descriptor.center)
            merge!(segments[seg.descriptor.center], Dict(seg.descriptor.target=>seg))
        else
            merge!(segments, Dict(seg.descriptor.center=>Dict(seg.descriptor.target=>seg)))
        end
    end

    g = NodeGraph{Int, Int}(SimpleGraph())
    ids = get_segments_naifids(segments)
    for id_i in unique(reduce(vcat, ids))
        add_vertex!(g, id_i)
    end 
    for idx in ids
        add_edge!(g, idx[1], idx[2])
    end
    SPK(daf, segments, g)
end

Base.show(io::IO, spk::SPK) = print(io, "SPK($(spk.segments[0][1].descriptor.name))")

segments(spk::SPK) = spk.segments

function list_segments(spk::SPK)
    s = String[]
    for (k,v) in spk.segments
        for l in keys(v)
            push!(s, "($k) => ($l)")
        end
    end
    sort!(s, lt=segstrlt)
end

function print_segments(spk::SPK)
    s = list_segments(spk)
    println(join(s, "\n"))
end

function segstrlt(a::String, b::String)
   rex = r"\([0-9]*\)$"
   ma = match(rex, a)
   mb = match(rex, b)
   ia = parse(Int, a[ma.offset+1:end-1])
   ib = parse(Int, b[mb.offset+1:end-1])
   ia < ib
end
