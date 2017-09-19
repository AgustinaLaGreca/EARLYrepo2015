function N = nstored(P);
% playlist/Nstored - number of samples stored in playlist object.
%
%   Nstored(P) or P.Nstored is the number of samples stored in playlist P, 
%   EXcluding the repititions. This is equal to sum(P.Nsam).
%
%   See also playlist, playllist/Nstored, playlist/SIZE, "methods playlist".

N = sum(Nsam(P));







