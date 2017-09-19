function P = not(P)
% playlist/not - set all reps of playlist to one.
%
%   not(P) or ~P is playlist P reduced by setting all reps to one.
%
%   See playlist/double, "methods playlist".

P.Nrep = 1+0*P.Nrep;





