function Bf = JLbflat(B);
% JLbflat - flatten binint output
for ii=1:numel(B),
    b = B(ii);
    bf.iexp = b.Yi.iGerbil;
    bf.icell = b.Yi.icell;
    bf.icond = b.Yi.icond;
    bf.freq  = b.Yi.Freq1;
    bf.SPL  = b.Yi.SPL1;
    bf.qqq_______________ = '____________';
    bf.Uidx_I = b.Yi.UniqueRecordingIndex;
    bf.Uidx_C = b.Yc.UniqueRecordingIndex;
    bf.Uidx_B = b.Yb.UniqueRecordingIndex;
    qq = rmfield(b ,{'Yi' 'Yc' 'Yb'});
    bf = structJoin(bf,qq);
    Bf(ii) = bf;
end