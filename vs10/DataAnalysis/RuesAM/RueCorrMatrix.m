function S = RueCorrMatrix(RS);
% RueCorrMatrix - grand matrix of cross-cell correlation
%   RueCorrMatrix(RS), RS is output from RueSelectSignCells
[MeanCorr, StdCorr, tt1p_right] = deal(nan(RS.Ncell));

Npair = numel(RS.CC);
for ipair=1:Npair,
    cc = RS.CC(ipair);
    idx1 = strmatch(cc.cell_1, RS.allNames);
    idx2 = strmatch(cc.cell_2, RS.allNames);
    [MeanCorr(idx1,idx2), MeanCorr(idx2,idx1)] = deal(cc.mean);
    [StdCorr(idx1,idx2), StdCorr(idx2,idx1)] = deal(cc.std);
    [tt1p_right(idx1,idx2), tt1p_right(idx2,idx1)] = deal(cc.tt1p_right);
end

DDiag = diag(MeanCorr);
disAtCorr = MeanCorr./kron(DDiag, DDiag.');

S = CollectInStruct(MeanCorr, StdCorr, disAtCorr, tt1p_right);


