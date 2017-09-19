function [Sn, Sb] = JLtestbin(ST, ExpID, iCell);
% JLtestbin - quick overview of binaurality of a given cell
%    JLtestbin(ST, ExpID, iCell);
%    Use
%      load D:\processed_data\JL\JLbeatStats\STall.mat
%    to load previous output ST of JLbeatstats('all').
if isnumeric(ExpID),
    switch ExpID
        case 78, ExpID = 'RG09178';
        case 91, ExpID =  'RG10191';
        case 95, ExpID =  'RG10195';
        case 96, ExpID =  'RG10196';
        case 97, ExpID =  'RG10197';
        case 98, ExpID =  'RG10198';
        case 01, ExpID =  'RG10201';
        case 04, ExpID =  'RG10204';
        case 09, ExpID =  'RG10209';
        case 14, ExpID =  'RG10214';
        case 16, ExpID =  'RG10216b';
        case 19, ExpID =  'RG10219';
        case 23, ExpID =  'RG10223';
        otherwise, error('Invalid ExpID');
    end
end

% select exp & cell
ii = strmatch(ExpID, {ST.ExpID}); 
Sn = ST(ii); 
Sn = Sn([Sn.icell]==iCell);

if isempty(Sn),
    error('No data having Nspike>3');
end
% binaural stuff: more than 3 spikes, VS stats
iB = strmatch('B', {Sn.StimType}); 
Sb = Sn(iB); 
Sb = Sb([Sb.Nspike]>3);
hist([Sb.VSbin],0:0.05:1)
xlabel('VSbin');
title([ExpID '/ ' num2str(iCell)]);


















