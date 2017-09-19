function [strSPL, strSpec1, strSpec2] = stimdispTHR(Stim);
% stimdispTHR - strings describing specific parameters of THR stimulus
%
%    See also stimGUI, dataset/stimlist, dataset/stimlist_strfun, databrowse.

STR = stimlist_strfun(dataset); % helpers for num->str conversion

strSPL = [STR.shstring(Stim.StartSPL) ':' STR.shstring(Stim.StepSPL) ':' STR.shstring(Stim.EndSPL) 'dB SPL']; 
strSpec1 = ''; 
strSpec2 = ''; 