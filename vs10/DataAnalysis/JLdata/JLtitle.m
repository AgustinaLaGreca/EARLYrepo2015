function TT = JLtitle(Uidx);
% JLtitle - plot title identifying a JL recording (single freq)
%    JLtitle(Uidx)

S = JLdbase(Uidx);
TT = sprintf('Exp %d cell %d  --- %s  %d Hz   %d dB  (series %d; icell_run=%d)', ...
    S.iexp, S.icell ,S.chan, S.freq, S.SPL, S.iseries, S.icell_run);



