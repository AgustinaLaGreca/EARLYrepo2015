function S = ctranspose(S);
% tsig/ctranspose - complex conjugate tsig object.
%    ctranspose(T) or T' is the complex conjugate of tsig T, obtained by
%    replacing all samples by their complex conjugate. Note that no
%    "transpose" operation can be performed on tsig objects, because they
%    do not have a matrix format.
%
%    See also tsig/ops, tsig/transpose.

S.Waveform = cellfun(@conj, S.Waveform, 'UniformOutput',false);



