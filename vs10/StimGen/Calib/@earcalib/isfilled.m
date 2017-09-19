function iv = isfilled(T);
% Earcalib/isfilled - true for a filled Earcalib object.
%    isvoid(T) returns true for a filled Earcalib object T, false otherwise.
%    A filled transfer object is one that has been measured by
%    Earcalib/measure.

iv = ~isempty(T.Transfer);






