function T = uminus(T);
% tsig/UMINUS - UMINUS for tsig objects
%   -T is T with all sample values negated.
%
%   See also tsig/uplus.

T.Waveform = cellfun(@uminus,T.Waveform, 'UniformOutput', false);



