function P = flip(P)
% playlist/flip - time-reverse playlist object
%
%   P = flip(P) reverses both the waveforms in P and their play order.
%
%   See playlist/not, "methods playlist".

P.Waveform = cellfun(@flipud, P.Waveform, 'UniformOutput', false);
P = list(P, fliplr(list(P)));



