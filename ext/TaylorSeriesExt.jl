module TaylorSeriesExt 

import Ephemerides: find_logical_record
using Ephemerides: DAF, 
                   SPKSegmentHeader2, SPKSegmentHeader8, SPKSegmentHeader20,
                   SPKSegmentHeader1, SPKSegmentHeader5, SPKSegmentHeader9,
                   SPKSegmentHeader14, SPKSegmentHeader18
using TaylorSeries: constant_term, Taylor1

function find_logical_record(head::SPKSegmentHeader2, time::Taylor1{<:Real})
    return find_logical_record(head, constant_term(time))
end

function find_logical_record(head::SPKSegmentHeader8, time::Taylor1{<:Real})
    return find_logical_record(head, constant_term(time))
end

function find_logical_record(head::SPKSegmentHeader20, time::Taylor1{<:Real})
    return find_logical_record(head, constant_term(time))
end

function find_logical_record(daf::DAF, head::SPKSegmentHeader1, time::Taylor1{<:Real})
    return find_logical_record(daf, head, constant_term(time))
end

function find_logical_record(daf::DAF, head::SPKSegmentHeader5, time::Taylor1{<:Real})
    return find_logical_record(daf, head, constant_term(time))
end

function find_logical_record(daf::DAF, head::SPKSegmentHeader9, time::Taylor1{<:Real})
    return find_logical_record(daf, head, constant_term(time))
end

function find_logical_record(daf::DAF, head::SPKSegmentHeader14, time::Taylor1{<:Real})
    return find_logical_record(daf, head, constant_term(time))
end

function find_logical_record(daf::DAF, head::SPKSegmentHeader18, time::Taylor1{<:Real})
    return find_logical_record(daf, head, constant_term(time))
end

end
