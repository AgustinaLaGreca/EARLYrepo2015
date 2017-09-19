function [strSPL, strSpec1, strSpec2] = stimdispZWSPL(Stim);
% stimdispZWSPL - strings describing specific parameters of ZWSPL stimulus
%    [strSPL, strSpec1, strSpec2] = stimdispZWSPL(Stim) returns strings
%    strSPL, strSpec1, strSpec2 describing the stimulus parameters of the
%    ZWSPL stimulus. Stim is the struct containing all the stimulus
%    parameters (see stimdefZWSPL). These strings are used by
%    dataset/stimlist and determine the listing in databrowse.
%
%    See also stimGUI, dataset/stimlist, dataset/stimlist_strfun, databrowse.

STR = stimlist_strfun(dataset); % helpers for num->str conversion

strSPL = '(per tone)';
strSpec1 = [STR.shstring([Stim.LowFreq Stim.HighFreq], '..') ' Hz']; % noise cutoffs
strSpec2 = [STR.shstring(Stim.ITD) ' ms ITD']; 


