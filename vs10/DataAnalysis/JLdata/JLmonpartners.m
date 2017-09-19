function [D3, P] = JLmonpartners(D);
% JLmonpartners - find monaural counterparts of binaural recordings
%    [D3, P] = JLmonpartners(D), where D is JLdbase output.
%       D3 is D restricted to binaural recordings with complete triple
%       P is struct describing partners of coresponding D3 elements

D3 = D([D.chan]=='B');
for ii=1:numel(D3),
    db = D3(ii);
    Dp = D([D.icell_run]==db.icell_run & [D.freq] == db.freq & [D.SPL] == db.SPL); % candidate partners
    dp.bin_partner = [Dp([Dp.chan]=='B').UniqueRecordingIndex];
    dp.ipsi_partner = [Dp([Dp.chan]=='I').UniqueRecordingIndex];
    dp.contra_partner = [Dp([Dp.chan]=='C').UniqueRecordingIndex];
    dp.multmon = numel(dp.ipsi_partner)>1 & numel(dp.contra_partner)>1;
    P(ii) = dp;
    bintrian(ii) = ~isempty(dp.ipsi_partner) && ~isempty(dp.contra_partner);
end
[D3, P] = deal(D3(bintrian), P(bintrian));


