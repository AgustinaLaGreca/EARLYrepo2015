function T = flipud(T)
% tsig/flipud - time-reverse tsig object
%
%   T = flipus(T) reverses the time-order of T.
%
%   See tsig, "methodshelp tsig", playlist/flip.

T.Waveform = cellfun(@flipud, T.Waveform, 'UniformOutput', false); % note that flipud time-reverses also any playlist waveforms


