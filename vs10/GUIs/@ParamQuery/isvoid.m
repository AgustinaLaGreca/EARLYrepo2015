function iv = disp(Q)
% ParamQuery/isvoid - true for void Paramquery objects.
%   isvoid(Q) returns when Q is a void ParamQuey. Arrays allowed.
%
%   See ParamQuery.

iv = cellfun('isempty', {Q.Name});
iv = reshape(iv, size(Q));