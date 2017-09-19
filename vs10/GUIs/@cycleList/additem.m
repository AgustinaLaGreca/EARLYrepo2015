function C = additem(C, Label, Userdata);
% CycleList/additem - add item to cyclelist
%   C=additem(C, Label, Userdata) adds a menuitem to Cycle list C.
%   The item will be placed on top of the list and has label Label.
%   When rendered by cyclelist/draw, the uimenu item's userdata property
%   will be set to Userdata, defaulting to [].
%   If an item having the same Label and Userdata is already present in the
%   list, additem places it on top of the list.  If C is already being
%   rendered, then its view is updated.
%
%   See also cycleList, cycleList/draw.

if nargin<3, Userdata=[]; end

if nargout<1, error('Too few output args'); end
newItem = CollectInStruct(Label,Userdata);
Nit = numel(C.Items);
% look if newItem already exists
ihit = [];
for ii=1:Nit, 
    if isequal(newItem, C.Items(ii)),
        ihit = ii;
        break;
    end
end
if ~isempty(ihit), % permute, don't add
    C.Items = C.Items([ihit 1:ihit-1 ihit+1:end]); % element is #ihit moved to the 1st place
else, % prepend new item; remove any item beyond the Nmax-th
    C.Items = [newItem C.Items];
    if numel(C.Items)>C.Nmax,
        C.Items = C.Items(1:C.Nmax);
    end
end

% refresh if rendered
if ishandle(C.Handle),
    refresh(C);
end







