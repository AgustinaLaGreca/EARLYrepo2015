function Y = duration(T)
% tsig/duration - duration time of tsig object
%   duration(T) or T.duration is the duration [ms] of T. For multi-channel
%   T, the duration is a 1 x T.nchan array.
%
%   Note: tsig/endtime is identical to tsig/offset.
%
%   See tsig/starttime, tsig.duration.

Y = dt(T).*nsam(T);
    


