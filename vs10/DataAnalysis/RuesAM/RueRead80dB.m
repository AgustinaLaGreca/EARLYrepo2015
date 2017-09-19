function [Y1, Y2, dt, Dur] = RueRead80dB(F1, F2, TimeWin);
% RueRead80dB(F1, F2, TimeWin)

persistent MasterList

if nargin<3, TimeWin = 450; end % 450 ms = whole trace

Ncond = 29; % first 29 conditions are the 80-dB ones, all freqs

if isempty(MasterList),
    [File, Nrep, Ncell] = ListRueData;
    MasterList = CollectInStruct(File, Nrep, Ncell);
end

[File, Nrep, Ncell] = deal(MasterList.File, MasterList.Nrep, MasterList.Ncell);
i1 = strmatch(F1, File);
i2 = strmatch(F2, File);

N1 = Nrep(i1);
N2 = Nrep(i2);

% read F1's
Y1 = 0;
for irep=0:N1-1,
    [Y, dt, Dur] = ReadRueData(F1, irep, Ncond, TimeWin);
    Y1 = Y1 + Y/N1;
end
% read F2's
Y2 = 0;
for irep=0:N2-1,
    [Y, dt, Dur] = ReadRueData(F2, irep, Ncond, TimeWin);
    Y2 = Y2 + Y/N2;
end












