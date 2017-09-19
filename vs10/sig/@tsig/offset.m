function Y = offset(T)
% tsig/offset - end time of tsig object
%   offset(T) or T.offset is the offset [ms] of T, which is a 1xT.nchan array. 
%   By definition, the offset equals the onset plus the duration.
%
%   Note: tsig/endtime is identical to tsig/offset.
%
%   See tsig/starttime, tsig.duration.

Y = T.t0+duration(T);
    


