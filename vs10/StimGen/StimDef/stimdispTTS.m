function [strSPL, strSpec1, strSpec2] = stimdispTTS(Stim);
% stimdispTTS - strings describing specific parameters of TTS stimulus
%    [strSPL, strSpec1, strSpec2] = stimdispTTS(Stim) returns strings
%    strSPL, strSpec1, strSpec2 describing the stimulus parameters of the
%    TTS stimulus. Stim is the struct containing all the stimulus
%    parameters (see stimdefTTS). These strings are used by
%    dataset/stimlist and determine the listing in databrowse.
%
%    See also stimGUI, dataset/stimlist, dataset/stimlist_strfun, databrowse.

STR = stimlist_strfun(dataset); % helpers for num->str conversion

strSPL = [STR.shstring(Stim.SPL) '-dB probe']; % probe level
strSpec1 = [STR.shstring(Stim.SupFreq) '-Hz sup'];  % suppr freq
strSpec2 = [STR.shstring(Stim.SupSPL) '-dB sup']; % suppr level





