function Mess=save(C);
% CycleList/save - save current state of CycleList object
%   Mess=save(C) attempts to save the current state of C.
%   By "state" is meant the list of Items present in C.
%   Mess is empty by default, but contains an error message in case of
%   problems. The filename depends on the name of the GUI in which C is 
%   placed (see CycleList/filename). Therefore, SAVE may only be applied to 
%   rendered CycleList objects. 
%
%   Typically, SAVE is evoked implicitly by REFRESH in order to store the
%   current state of C.
%
%   See also cycleList, cycleList/filename, cycleList/load.

figh = parentfigh(C.Handle);
figh = figh.Number;
if ~ishandle(figh),
    error('CycleList is not rendered - cannot save.');
end
GUIname = getGUIdata(figh,'GUIname');
FN = filename(C,GUIname);

Citems = C.Items;
try,
    save(FN, 'Citems', '-mat');
    Mess = '';
catch,
    Mess = lasterr;
end






