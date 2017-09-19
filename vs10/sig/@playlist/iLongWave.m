function I = iLongWave(P,k);
% playlist/iLongWave - expanded waveform index vector of playlist object.
%
%   iLongWaveg(P) or P.iLongWave returns the expanded version of the waveform 
%   indices to be played, meaning that the repetitions are "written out."
%
%   iLongWave(P,K) only returns the Kth index or indices.
%
%   See playlist, playlist/iWave, Playlist/list.

I = [];
for ii=1:length(P.iplay),
    iwa = P.iplay(ii);
    nrep = P.Nrep(ii);
    I = [I repmat(iwa,1,nrep)];
end

if nargin>1,
    I = I(k);
end




