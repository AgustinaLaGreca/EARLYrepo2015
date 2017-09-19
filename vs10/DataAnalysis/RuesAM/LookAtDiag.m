function Diag = LookAtDiag(Twin);
% LookAtDiag - select diagional correlations
FN = ['AllStats_' num2str(Twin) 'ms.mat'];
FFN = fullfile(['D:\Data\RueData\IBW\Audiograms\Corr_' num2str(Twin)  'ms'], FN);
qq = load(FFN);
qq = qq.Stats;
nam1 = {qq.cell_1};
nam2 = {qq.cell_2};
nam=unique(nam1);
Ncell = numel(nam);
idiag = [];
for icell=1:Ncell,
    cc = nam{icell};
    ihit1 = strmatch(cc, nam1, 'exact');
    ihit2 = strmatch(cc, nam2, 'exact');
    idiag = [idiag , intersect(ihit1, ihit2)];
end
Diag = qq(idiag);








