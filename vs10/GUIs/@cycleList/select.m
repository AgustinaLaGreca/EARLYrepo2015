function Select(src, dum, void, CLhandle, i_item); 
% CycleList/Select - generic callback of cycle list
%   Select is the generic callback function of CycleList menu items. 
%   Select delegates the work by passing four arguments: Src, Event, L and 
%   Item to the callback function of the clicked cycleList menu item:
%       feval(CB, Src, Event, L, Item)
%   CB is the callback function of L.
%   Src is the handle of the the calling menu item.
%   Event is some future Matlab thingy.
%   L is the current state of the CycleList object whose item was selected.
%   Item is the selected item from the array L.Items.
%
%   An sample implementation is found in StimGuiFilePulldown.
%
%   See CycleList, StimGuiFilePulldown.

dumdum={src void dum}; % to silence Mlint ;)

% retrieve calling cycle list
C = get(CLhandle, 'userdata');

% which item was calling?
item = C.Items(i_item);

% call true callback with the GUI's handle & item's userdata as arguments
feval(C.Callback, src, dum, C, item);

% % place item on top of stack
% C = additem(C,item.Label, item.Userdata);
