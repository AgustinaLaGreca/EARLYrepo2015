function iv = isvoid(P);
% earcalib/isvoid - true for a void earcalib object.
%   earcalib(PTC) returns True if PTC is a void earcalib, false
%   otherwise.

iv = isequal([], P.Experiment);