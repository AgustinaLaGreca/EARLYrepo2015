function N = nsamplay(P)
% playlist/Nsamplay - total sample count of playlist object.
%
%   Nsamplay(P) or P.Nsamplay is the total number of samples to be played as
%   when playing playlist P, including the repititions. 
%   The number of samples stored in P is generally lower, because multiple
%   reps are stored only once.
%
%   See also playlist, playllist/Nstored, playlist/SIZE, "methods playlist".

Nsam = nsam(P); % # samples in stored waveforms
N = sum(nrep(P).*Nsam(P.iplay));







