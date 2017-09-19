function Y = starttime(T, t0)
% tsig/starttime - get/set onset of tsig object
%
%   starttime(T) or T.starttime is the onset [ms] of T, which is a 1xT.nchan array. 
%
%   T=starttime(T,t0) or T.starttime=t0 sets the onset to t0 ms.
%   t0 is either a scalar or a row 1xT.nchan array, whose elements set the
%   onsets of the individual channels. Note that setting the starttime of
%   a single-channel T to an 1xM array value transforms T into an M-channel
%   tsig.
%
%   Note. tsig/starttime is identical to tsig/onset.
%
%   See tsig/onset, tsig/endtime, tsig.duration.

if nargin==1, % get
    Y = T.t0;
elseif nargin==2, % set
    T.t0 = t0;
    [T, Mess] = test(T); % let tsig constructor do the checking
    error(Mess);
    Y = T;
else, 
    error('Invalid #input arguments.');
end
    


