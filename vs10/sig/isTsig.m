function IP = isTsig(P);
% isTsig - true for tsig object.
%    isTsig(T) returns logical true (1) if T is a tsig object and
%    logical false (0) otherwise.
%
%    See also tsig.

IP = isa(P, 'tsig');