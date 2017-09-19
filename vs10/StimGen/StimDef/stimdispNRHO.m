function [strSPL, strSpec1, strSpec2] = stimdispNRHO(Stim);
% stimdispNRHO - strings describing specific parameters of NRHO stimulus
%    [strSPL, strSpec1, strSpec2] = stimdispNRHO(Stim) returns strings
%    strSPL, strSpec1, strSpec2 describing the stimulus parameters of the
%    NRHO stimulus. Stim is the struct containing all the stimulus
%    parameters (see stimdefNRHO). These strings are used by
%    dataset/stimlist and determine the listing in databrowse.
%
%    See also stimGUI, dataset/stimlist, dataset/stimlist_strfun, databrowse.

STR = stimlist_strfun(dataset); % helpers for num->str conversion

strSPL = [STR.shstring(Stim.SPL) ' ' Stim.SPLUnit];
strSpec1 = STR.modstr(Stim); 
strSpec2 = [STR.shstring([Stim.LowFreq Stim.HighFreq], '..') ' Hz']; % noise cutoffs






