function [MC, Stim, CMX] = JSmcoutput(ExpID, iRec);
% JSmcoutput - cached maskcorr output
%   EXAMPLE
%   [MC, Stim, CondMetrics] = JSmcoutput('RG11295',34); output of maskcorr('RG11295',34)
%   [MC, Stim, CondMetrics] = JSmcoutput('RG11295',34:36); cell array of outputs

if numel(iRec)>1, % recursion
    MC = {}; Stim = []; CMX = [];
    for ii=1:numel(iRec),
        try, [mc, st, cmx] = JSmcoutput(ExpID, iRec(ii));
        catch, continue;
        end
        MC = [MC, mc];
        Stim = [Stim st];
        CMX = [CMX, cmx];
    end
    return;
end

%===single iRec from here=============

qq = dir(which('dataset/maskcorr')); % retrieve date of maskcorr.m
[MC_Stim, CFN, CP] = getcache(mfilename, {ExpID, iRec qq.datenum}); % new date -> new computation
if ~isempty(MC_Stim),
    MC = MC_Stim.MC;
    CMX = MC_Stim.CMX;
    Stim = MC_Stim.Stim;
    return;
end

ds = read(dataset(), ExpID, iRec);
[MC, CMX] = maskcorr(ds);
delete(gcf);
Stim  = ds.stim;
putcache(CFN, 5000, CP, CollectInStruct(MC, Stim, CMX));




