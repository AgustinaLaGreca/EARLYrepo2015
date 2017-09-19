function [FN, Ex] = filename(C, GUIname);
% CycleList/filename - full filename for storage and retrieval of CycleList
%   [FN, Exists] = filename(C, GUIname) returns the full filename FN of the
%   file in which to store CycleList C when placed in a GUI named GUIname.
%   The output argument Exists is a logical indicating whether FN exists.
%   Filename is used by the LOAD and SAVE methods for CycleList objects.
%
%   See also cycleList, cycleList/load, cycleList/save, GUIdefaultsDir.

DD = GUIdefaultsDir(GUIname);

FN = fullfile(DD,[GUIname '_' C.Name '.CycleList']);
Ex = exist(FN,'file');






