function P = PulldownMenuContaining(C)
% CycleList/PulldownMenuContaining - PulldownMenu Containing cycleList
%   PulldownMenuContaining(C) returns the PulldownMenu object in which C is
%   contained. C must be rendered.
%
%   See Pulldown/additem, CycleList, PulldownMenu/getCycleList.

hP = C.parentHandle;
if ~isSingleHandle(hP),
    error('Cannot find containing PulldownMenu of non-rendered CycleList');
end
P = get(hP,'userdata');




