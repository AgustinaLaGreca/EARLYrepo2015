function [S, ST] = DSinfoString(ds);
% SGSRdataset/DSinfoString - summary of ID and parameters of dataset
%   S = DSinfoString(DS) is a 1x7 cell array of strings
%   summarizing the ID and parameters of dataset DS.
%   All strings contain a trailing space to facilitate GUI rendering.
%
%   The following table explains what the strings S{k} contain:
%     k       content           Example     comment
%    ------------------------------------
%     1   dataset seq n          -42     (negative means IDF/SPK)
%     2   seq ID string       <13-22-FS>  (<stimtype> if no seq ID)
%     3     carrier freq        200/250 Hz
%         or noise cutoffs    100-3000 Hz
%     4        SPL              30/40 dB
%     5     stimulus-specifics   F, rho=0.95  (noise params)
%                              25 Hz, 100%  (modulation params)
%                                100 us       (click duration)
%                              myWavList   (filename of wavlist)
%     6    timing      4* x 20 x 1000/1300 ms (asterisk: incomplete rec)
%     7    DAchan                B
%
%   [S, ST] = DSinfoString(DS) also returns a struct ST whose fields are
%   the cells os S, viz: iseqStr, IDstr, Fstr, SPLstr, SPstr, Tstr, DAstr.
%
%   See also showStimParam.

% delegate to non-method fcn that is allowed to use SGSRdataset/subref
[S, ST] = nonmethod_DSinfoString(ds);
