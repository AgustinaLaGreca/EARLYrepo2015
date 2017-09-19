function Y = onset(T, t0)
% tsig/onset - get/set onset of tsig object
%
%   onset(T) or T.onset is the onset [ms] of T, which is a 1xT.nchan array. 
%
%   T=onset(T,t0) or T.onset=t0 sets the onset to t0 ms.
%   t0 is either a scalar or a row 1xT.nchan array, whose elements set the
%   onsets of the individual channels. Note that setting the onset of
%   a single-channel T to an 1xM array value transforms T into an M-channel
%   tsig.
%
%   Note. tsig/onset is identical to tsig/starttime.
%
%   See tsig, "methodshelp tsig".

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
    


