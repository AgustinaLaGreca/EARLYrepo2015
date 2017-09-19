function icr = JL_icell_run(Uidx)
% JL_icell_run - icell_run index of recording
%    JL_icell_run(Uidx)

persistent LUT
if isempty(LUT),
    DB = JLdbase;
    LUT = [[DB.UniqueRecordingIndex]; [DB.icell_run]];
end

icr = LUT(2, ismember(LUT(1,:),Uidx));




