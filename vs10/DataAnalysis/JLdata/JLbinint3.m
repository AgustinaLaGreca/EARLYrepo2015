function B = JLbinint3(Utix, doPlot);
% JLbinint3 - analyze bin. interaction from I/C/B data
%   JLbinint3(Utix), where Utix is unique triplet index, wich equals minus
%   the Uidx of the B recording.
%
%   See also JLbinint2.

SV = retrieve_dbase('D:\processed_data\JL\JLsplitvar\JLsplitvar.dbase', Utix);

BI= JLbinint(JLdatastruct(SV.Uidx_I), JLdatastruct(SV.Uidx_C), JLdatastruct(SV.Uidx_B));


B = SV;



1+1
















