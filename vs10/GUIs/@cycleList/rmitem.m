function C = rmitem(C, k);
% CycleList/rmitem - remove item(s) from cyclelist
%   C=rmitem(C, k) removes the k-th item from C. When k is an array, 
%   multiple items are removed. If C is already being
%   rendered, then its view is updated.
%   
%   C=rmitem(C, 'Foo') removes the item whose label is Foo. Exact,
%   case-sensitive match required.
%   
%   C=rmitem(C) removes all items of C.
%
%   See also cycleList, cycleList/draw.

if nargin<2, k=1:length(C.Items); end


if ischar(k),
    k = find(strcmp(k,{C.Items.Label}));
end
C.Items(k) = [];

% refresh if rendered
if ishandle(C.Handle),
    refresh(C);
elseif nargout<1, 
    error('Too few output args; cycle list not rendered.');
end





