function m=RuePlotShiftCorr80dB(f1, f2, TimeWin);
% Usage: RuePlotShiftCorr80dB(f1, f2, TimeWin);

if nargin<3, TimeWin=450; end

global M
if isempty(M),
    FN = ['D:\Data\RueData\IBW\Audiograms\RueShiftCorr80dB_' num2str(TimeWin) 'ms.mat'];
    M=load(FN);
    M = M.M(:);
end

i1 = strmatch(lower(f1), lower({M.f1}));
i2 = strmatch(lower(f2), lower({M.f2}));
ihit = intersect(i1,i2);
if isempty(ihit),
    i1 = strmatch(lower(f1), lower({M.f2}));
    i2 = strmatch(lower(f2), lower({M.f1}));
    ihit = intersect(i1,i2);
end
m = M(ihit(1));

[Y1, Y2, dt, Dur] = RueRead80dB(m.f1, m.f2, TimeWin);
NsamCond = round(Dur/29/dt);
localPlot(m, dt, Y1, Y2, NsamCond, m.NcondOverlap);

function [y1, y2] = localPlot(m, dt, Y1, Y2, NsamCond, NcondOverlap);
Ncond = round(numel(Y1)/NsamCond);
MaxShift = Ncond-NcondOverlap;
% MAX
[y1, y2] = localSelect(Y1,Y2, m.OffsetMax(1), m.OffsetMax(2), NcondOverlap, NsamCond);
subplot(2,1,1);
dplot(dt, y1, 'color', [0 0 0.7]);
xdplot(dt, y2, 'color', [0 0.7 0]);
title([m.f1 ' vs. ' m.f2 ' C_{max}=' num2str(m.Cmax,3)]);
ylabel('V_m (mV)');
cond1 = m.OffsetMax(1)+[1 NcondOverlap];
cond2 = m.OffsetMax(2)+[1 NcondOverlap];
legend({['cond ' num2str(cond1(1)) '..' num2str(cond1(2))], ['cond ' num2str(cond2(1)) '..' num2str(cond2(2))]}, ...
    'location', 'northwest');
% MIN
[y1, y2] = localSelect(Y1,Y2, m.OffsetMin(1), m.OffsetMin(2), NcondOverlap, NsamCond);
subplot(2,1,2);
dplot(dt, y1, 'color', [0 0 0.7]);
xdplot(dt, y2, 'color', [0 0.7 0]);
title([' C_{min}=' num2str(m.Cmin,3)]);
ylabel('V_m (mV)');
xlabel('Time (ms)');
cond1 = m.OffsetMin(1)+[1 NcondOverlap];
cond2 = m.OffsetMin(2)+[1 NcondOverlap];
legend({['cond ' num2str(cond1(1)) '..' num2str(cond1(2))], ['cond ' num2str(cond2(1)) '..' num2str(cond2(2))]}, ...
    'location', 'northwest');

function [y1, y2] = localSelect(Y1,Y2,ioffset1, ioffset2, NcondOverlap, NsamCond);
isam1_start = 1 + ioffset1*NsamCond;
isam1_end = isam1_start + NcondOverlap*NsamCond - 1;
y1 = Y1(isam1_start:isam1_end);

isam2_start = 1 + ioffset2*NsamCond;
isam2_end = isam2_start + NcondOverlap*NsamCond - 1;
y2 = Y2(isam2_start:isam2_end);














