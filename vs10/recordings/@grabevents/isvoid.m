function v=isvoid(E);
% Eventdata/isvoid - true for void Eventdata object.
%   isvoid(E) returns TRUE if E is a void Eventdata object, FALSE otherwise
%   A void Eventdata object is the object returned by the Eventdata 
%   constructor called without any input arguments. 
%
%   See Eventdata.

for ii=1:numel(E),
    v(ii) = isempty(E(ii).Fsam);
end
v = reshape(v, size(E));
