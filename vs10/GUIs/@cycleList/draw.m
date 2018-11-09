function C=draw(h, C, XY);
% CycleList/draw - draw cycleList object.
%   draw(h,C) renders cycleList C in the uimenu control with handle h.
%
%   See CycleList, CycleList/refresh, GUIpiece/draw, GUIpiece.

C.parentHandle = h;
% create invisible uimenu that serves as "container"
C.Handle = uimenu(C.parentHandle,...
        'visible', 'off', ...
        'tag', ['Cycle list ' C.Name]);

% Create uimenuitems
Sepa = 'on';
for ii=1:C.Nmax
    C.ItemHandles(ii) = uimenu(C.parentHandle, ...
        'visible', 'off', ...
        'separator', Sepa, ...
        'Callback', {@cyclelist_callback, C, ii}, ... % the first arg serves to invoke the select method
        'Label', ['unused item #' num2str(ii)], ...
        C.uimenuProps);
    Sepa = 'off'; % only first item has separator
end

% Store C in GUIdata of figure
figh = parentfigh(h);
figh = figh.Number;
CL = getGUIdata(figh, 'CycleList', []);
CL = [CL C];
setGUIdata(figh, 'CycleList', CL);

% if no items are present, try to retrieve them from file
if isempty(C.Items), 
    C = load(C);
end

% now synchronize the rendering with the items in C.Items. This is
% delegated to the refresh method, which also updates the GUIdata of the
% figure and the userdata of the "container" iumenu with handle C.Handle.
C=refresh(C);

% function that replaces @cyclelist/select
function cyclelist_callback(src, event, C, i_item)
% CycleList_callback - generic callback of cycle list
%   Generic callback function of CycleList menu items. 
%   Select delegates the work by passing four arguments: Src, Event, C and 
%   Item to the callback function of the clicked cycleList menu item:
%       feval(CB, Src, Event, L, Item)
%
%   Src is the handle of the calling menu item.
%   Event is necessary for callbacks definitions.
%   C is the initial cyclelist object whose item was selected, empty (no saved files)
%   CB is the callback function of CL.
%   CL is the actual up to date cyclelist object, obtained from C as
%   get(C.Handle,'userdata')
%   Item is the selected item from the array L.Items.
%
%   An sample implementation is found in StimGuiFilePulldown.
%
%   See CycleList, StimGuiFilePulldown.

% retrieve calling cycle list current values
C = get(C.Handle, 'userdata');

% which item was calling?
item = C.Items(i_item);

% call true callback with the GUI's handle & item's userdata as arguments
feval(C.Callback, src, event, C, item);
