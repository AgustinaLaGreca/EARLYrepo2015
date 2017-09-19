function v=isvoid(G);
% grab_adc/isvoid - true for void Grab_adc object.
%   isvoid(G) returns TRUE if E is a void Grab_adc object, FALSE otherwise
%   A void Grab_adc object is the object returned by the Grab_adc 
%   constructor called without any input arguments. 
%
%   See Eventdata.

for ii=1:numel(G),
    v(ii) = isempty(G(ii).Fsam);
end
v = reshape(v, size(G));
