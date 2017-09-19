function A = all(S);
% tsig/all - ALL for tsig objects
%    all(S) performs a channel-wise ALL on the samples of tsig object S and
%    returns the logical result in a 1 x S.nchan  row vector.
%
%    See also tsig/any, ALL.

A = cellfun(@all, S.Waveform);



