function A = any(S);
% tsig/any - ANY for tsig objects
%    any(S) performs a channel-wise ANY on the samples of tsig object S and
%    returns the logical result in a 1 x S.nchan  row vector.
%
%    See also tsig/all, ANY.

A = cellfun(@any, S.Waveform);



