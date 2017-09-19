function M = EvalRueShiftCorr80dB(TimeWin);
% Usage: function EvalRueShiftCorr80dB(TimeWin);

if nargin<1, TimeWin=450; end

FN = ['D:\Data\RueData\IBW\Audiograms\RueShiftCorr80dB_' num2str(TimeWin) 'ms.mat'];

M=load(FN);
M = M.M;  % matrix -> col vector
Nrec = round(sqrt(numel(M)));

% fix flipped storage of offsets when i1<i2
for i2=2:Nrec,
    for i1=1:i2-1,
        M(i1,i2).OffsetMin = M(i1,i2).OffsetMin([2 1]);
        M(i1,i2).OffsetMax = M(i1,i2).OffsetMax([2 1]);
    end
end

for i1=1:Nrec,
    for i2=1:Nrec,
        m = M(i1,i2);
        OffsetMin1(i1,i2) = m.OffsetMin(1);
        OffsetMin2(i1,i2) = m.OffsetMin(2);
        OffsetMax1(i1,i2) = m.OffsetMax(1);
        OffsetMax2(i1,i2) = m.OffsetMax(2);
    end
end

M = M(:);  % matrix -> col vector
SqSize = Nrec*[1 1];
Cmin = reshape([M.Cmin], SqSize);
Cmax = reshape([M.Cmax], SqSize);

% ====== save ==========
Prefix = ['D:\Data\RueData\IBW\Audiograms\RueShiftCorr80dB_' num2str(TimeWin) 'ms_'];
FN1 = [Prefix, 'Cmax.xls']; xlswrite(FN1, Cmax);
FN2 = [Prefix, 'Cmin.xls']; xlswrite(FN2, Cmin);
FN3 = [Prefix, 'OffsetMax1.xls']; xlswrite(FN3, OffsetMax1);
FN4 = [Prefix, 'OffsetMax2.xls']; xlswrite(FN4, OffsetMax2);
FN5 = [Prefix, 'OffsetMin1.xls']; xlswrite(FN5, OffsetMin1);
FN6 = [Prefix, 'OffsetMin2.xls']; xlswrite(FN6, OffsetMin2);

RecNames  = {M(1:Nrec).f1}.';
xlswrite([Prefix 'CellNames.xls'], RecNames);








