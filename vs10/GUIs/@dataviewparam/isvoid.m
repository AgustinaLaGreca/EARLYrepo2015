function iv = isvoid(P);
% dataviewparam/isvoid - true for void dataviewparamobjacte.
%   isvoid(P) returns True when P is void; False otherwise. A void P is a P
%   with [] as Dataviewer.

iv = isequal(@nope, P.Dataviewer);




