function dt = dt(T)
% tsig/dt - sample period of tsig object
%
%   dt(T) or T.dt is the sample period [ms] of T.
%
%   See tsig, "methodshelp tsig".

dt = 1/T.Fsam; % kHz -> ms


