function [strSPL, strSpec1, strSpec2] = stimdispITDF(Stim);
% stimdispITDF - strings describing specific parameters of ITDF stimulus
%    [strSPL, strSpec1, strSpec2] = stimdispITDF(Stim) returns strings
%    strSPL, strSpec1, strSpec2 describing the stimulus parameters of the
%    ITDF stimulus. Stim is the struct containing all the stimulus
%    parameters (see stimdefITDF). These strings are used by
%    dataset/stimlist and determine the listing in databrowse.
%
%    See also stimGUI, dataset/stimlist, dataset/stimlist_strfun, databrowse.

STR = stimlist_strfun(dataset); % helpers for num->str conversion

strSPL = [STR.shstring(Stim.SPL) ' ' Stim.SPLUnit];
strSpec1 = [STR.shstring(Stim.startITD) ':' STR.shstring(Stim.stepITD) ':' STR.shstring(Stim.endITD) ' ms']; % ITD range
strSpec2 = STR.modstr(Stim); 






