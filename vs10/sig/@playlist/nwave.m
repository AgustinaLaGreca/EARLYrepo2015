function N = nwave(P)
% playlist/nwav - number of waveforms contained in playlist object.
%
%   nwave(P) or P.nwave returns the number of waveforms stored in playlist item P.
%
%   See playlist, playlist/nsam.

N = length(P.Waveform);






