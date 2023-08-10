using SPICE
using Test 

include("../spk/reader.jl")
include("../spk/spk1.jl")

filename = "res/example1spk_seg1.bsp";   

array = mmap(filename, Vector{UInt8});
header = load_daf(array)

summaries = get_daf_summaries(array, header);
for summ in summaries 
    println(parse_spk_segment_header(summ, header.lend))
end

spkhead = parse_spk_segment_header(summaries[1], header.lend)
segment = SPKSegmentType1(array, spkhead, header.lend)


find_logical_record(array, spkhead, segment, true, 0.0);

function find_logical_record(array, head::SPKSegmentHeader, seg::SPKSegmentType1, 
            lend::Bool, time::Number)

    i0 = 8*(head.faa-1)

    # search target epoch in epoch directories
    idx, found = 0, false 
    while !found && idx <= seg.ndirs 
        idx += 1

        epoch_dir = get_float(array, i0 - 8*(seg.ndirs - idx + 1), lend)
        if epoch_dir >= time 
            found = true 
        end

    end

    # Compute epoch table start and stop indexes! 
    if found 
        start_idx = (idx - 1) * 100 + 1
        stop_idx = 100*idx
    else
        start_idx = 100*seg.ndirs + 1
        stop_idx = seg.n
    end

    # GET LOGICAL RECORD WITHIN THAT DIRECTORY

    # i0 per epoch table (dopo segmenti record)
    i0 = 8*(head.iaa - 1) + 8*71*seg.n

    # i1 = i0 + 8*(start_idx - 1)

    record_idx, found = start_idx-1, false
    while record_idx <= stop_idx && !found 
        record_idx += 1 
        epoch_time = get_float(array, i0 + 8*(record_idx-1), lend)
        println(epoch_time)
        if epoch_time >= time 
            found = true 
        end
    end

    return found ? record_idx : stop_idx 

    # upper_boundary = get_float(array, i0 + 8*(record_idx-1), header.lend)
    # lower_boundary = record_idx != 1 ? 
                # get_float(array, i0 + 8*(record_idx-2), header.lend) : 
                # spkhead.tstart


end


# TESTING 
furnsh("/home/michele/spice/kernels/lsk/naif0012.tls", filename)

for i = 1:5000 
    et = rand(spkhead.tstart:spkhead.tend)
    et = spkhead.tend

    p3 = spkpos(string(spkhead.tid), et, "J2000", "NONE", string(spkhead.cid))[1]
    pos = vector3(array, header, spkhead, segment, et)

    @test p3 ≈ pos atol=1e-14 rtol=1e14

end

