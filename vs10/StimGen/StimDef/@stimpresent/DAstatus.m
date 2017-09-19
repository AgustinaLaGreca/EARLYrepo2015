function S=DAstatus(SP,isam)
% stimPresent/DAstatus - status of ongoing stimulus presentation
%    DAstatus(SP,isam), where isam is the sample index of the sequenced D/A
%    conversion, and SP is stimpresent object used for the stimulus 
%    presentation, returns a struct containing the following fields:
%         icond: index of stimulus condition now being played
%          irep: rep count of that stimulus condition 
%       remtime: remaining time [ms] of entire presentation
%
%   If isam is not specified, the 'isam_abs' value returned by 
%   SeqplayStatus is used.
%
%   See also stimpresent, sortConditions, seqplaystatus, DACprogress.

if nargin<2,
    Stat = seqplaystatus;
    isam = Stat.isam_abs;
end

[dum, S.ipres] = NthFloor(isam, SP.SamOffset);
if isequal(0, S.ipres), % before stim onset
    S.icond=-inf; 
    S.irep = 0;
    S.remtime = SP.TotDur;
elseif isequal(SP.Npres+1, S.ipres), % after stim offset
    S.icond=inf;
    S.irep = 0;
    S.remtime = 0;
else, 
    S.icond = SP.iCond(S.ipres,1);
    S.irep = SP.iRep(S.ipres);
    S.remtime = SP.TotDur-1e3*isam/SP.Fsam;
end






