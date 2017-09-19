function it = istoggle(Q)
% ParamQuery/istoggle - true for toggle Paramquery objects.
%   istoggle(Q) returns when Q is a toggle ParamQuey. Arrays allowed.
%   A paramquery is considered a toggle if its Unit property is a cellstring
%   and its String property is empty.
%
%   See ParamQuery.

it = cellfun(@iscellstr, {Q.Unit}) & cellfun(@isempty, {Q.String});
it = reshape(it, size(Q));



