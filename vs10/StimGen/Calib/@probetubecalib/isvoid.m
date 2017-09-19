function iv = isvoid(P);
% Probetubecalib/isvoid - true for a void Probetubecalib object.
%   Probetubecalib(PTC) returns True if PTC is a void Probetubecalib, false
%   otherwise.

iv = isequal([], P.Tube);