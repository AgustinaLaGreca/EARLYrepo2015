function M = RueShiftCorr80dB(NcondOverlap, TimeWin);
% Usage: M = RueShiftCorr80dB(NcondOverlap, TimeWin);

if nargin<1, NcondOverlap = 4; end
if nargin<2, TimeWin=450; end
[File, Nrep, Ncell] = ListRueData;
more off
Nrec = numel(File); % # recordings, including ipsi/contra

FN = ['D:\Data\RueData\IBW\Audiograms\RueShiftCorr80dB_' num2str(TimeWin) 'ms.mat'];
for i1=1:Nrec,
    for i2=1:i1,
        f1 = File{i1}, f2 = File{i2}
        [Y1, Y2, dt, Dur] = RueRead80dB(f1, f2, TimeWin);
        NsamCond = round(TimeWin/dt);
        [Cmin, Cmax, OffsetMin, OffsetMax] = localBestCorr(Y1, Y2, NsamCond, NcondOverlap);
        M(i1,i2) = CollectInStruct(f1, f2, Cmin, Cmax, OffsetMin, OffsetMax, NcondOverlap, TimeWin);
        M(i2,i1) = CollectInStruct(f1, f2, Cmin, Cmax, OffsetMin, OffsetMax, NcondOverlap, TimeWin);
        save(FN, 'M');
    end
end

% ====================================================
function [Cmin, Cmax, OffsetMin, OffsetMax] = localBestCorr(Y1, Y2, NsamCond, NcondOverlap);
Ncond = round(numel(Y1)/NsamCond);
MaxShift = Ncond-NcondOverlap;
allCorr = [];
OFF =[];
for ioffset1 = 0:MaxShift,
    isam1_start = 1 + ioffset1*NsamCond;
    isam1_end = isam1_start + NcondOverlap*NsamCond - 1;
    y1 = Y1(isam1_start:isam1_end);
    for ioffset2 = 0:MaxShift,
        isam2_start = 1 + ioffset2*NsamCond;
        isam2_end = isam2_start + NcondOverlap*NsamCond - 1;
        OFF = [OFF, [ioffset1; ioffset2]];
        y2 = Y2(isam2_start:isam2_end);
        allCorr = [allCorr, corr(y1,y2)];
    end
end
[Cmin, imin] = min(allCorr);
[Cmax, imax] = max(allCorr);
OffsetMin = OFF(:, imin);
OffsetMax = OFF(:, imax);








