function [Y, Nrep, dt] = readRueCond(FN, iSPL, iFreq);
% readRueCond - read all reps of Rue data, single condition
%   [Y, Nrep, dt] = readRueCond(FN, iSPL, iFreq);
%         Y: matrix whose columns are the reps
%      Nrep: # reps
%        dt: sample period in ms
%
%   See also ReadRueRep.

Nrep = RueNrep(FN);
NsamCond = 11250; % # samples per condition
Nfreq = 29;
icond = iFreq + Nfreq*(iSPL-1);
Y = [];
for irep=1:Nrep,
    [y, Ncond, dt] = readRueRep(FN, irep);
    Y = [Y, y(:, icond)];
end
