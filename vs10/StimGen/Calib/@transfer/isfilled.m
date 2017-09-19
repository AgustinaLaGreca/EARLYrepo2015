function iv = isfilled(T);
% transfer/isfilled - true for a filled Transfer object.
%    isvoid(T) returns true for a filled Transfer object T, false otherwise.
%    A filled transfer object is one that has been measured by
%    Transfer/measure.

iv = ~isempty(T.Ztrf);






