function P = real(P)
% playlist/real - REAL for playlist objects.
%
%   See playlist.

P.Waveform = cellfun(@real, P.Waveform, 'UniformOutput', false);







