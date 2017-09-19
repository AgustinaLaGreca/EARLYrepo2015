function S = transpose(S);
% tsig/transpose - no transpose for tsig objects.
%    transpose(T) or T.' is an error, because tsig objects do not have a 
%    matrix format. T' is allowed, and is the same as conj(T).
%
%    See also tsig/conj, tsig.transpose.

S=[]; % to please Mlint ;)
error('Transpose of tsig object is invalid.')


