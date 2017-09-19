function Select(src,dum, dum2, Cname, i_item)
% CycleList/Select - generic callback of cycle list
%   Extent(CL) or CL.Extent retirns [0 0].
%
%   See CycleList.

%src,dum, dum2, Cname, i_item, clickSide

% retrieve calling cycle list
figh = parentfigh(src);
C = getGUIdata(figh, 'CycleList');
iC = find(C,Cname); C= C(iC);
struct(C)
% which item was calling?
item = C.Items(i_item)

% call true callback with the GUI's handle & item's userdata as arguments
feval(C.Callback, figh, item.Userdata);

% place item on top of stack
C = additem(C,item.Label, item.Userdata);
