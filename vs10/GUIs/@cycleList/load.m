function C = load(C);
% CycleList/load - retrieve CycleList items from file.
%   C=load(C) attempts to load the previously saved state of C.
%   By "state" is meant the list of Items present when C was last saved.
%   The filename depends on the name of the GUI in which C is placed
%   (see CycleList/filename). Therefore, LOAD may only be applied to 
%   rendered CycleList objects. 
%
%   Typically, LOAD is evoked implicitly by draw in order to restore the
%   previous state of C.
%
%   See also cycleList, cycleList/filename, cycleList/save.

if nargout<1, error('Too few output args'); end
figh = parentfigh(C.Handle);
figh = figh.Number;
if ~ishandle(figh),
    error('CycleList is not rendered - cannot load.');
end
GUIname = getGUIdata(figh,'GUIname');
[FN, Exists] = filename(C,GUIname);

if Exists, 
    qq=load(FN, '-mat');
    C.Items = qq.Citems;
end




