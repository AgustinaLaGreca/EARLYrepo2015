function J = JL_jercog(Uidx);
% JL_jercog - check Jercog hypothesis on EPSP slope asymmetry
%    database name: 
%    DBFN = fullfile(processed_datadir, 'JL', 'Jercog', 'JL_jercog.dbase')
%       (on SIUT: D:\processed_data\JL\Jercog\JL_jercog.dbase)

if nargin<1, % all binaural recs
    DB = structMerge('UniqueRecordingIndex', JLdbase, JL_rec_select);
    DB = DB([DB.CanBeUsed]);
    DB = DB([DB.chan]=='B');
    Uidx = [DB.UniqueRecordingIndex];
end

% multiple recs: recursive call
if numel(Uidx)>1,
    for ii=1:numel(Uidx),
        ii
        J(ii) = JL_jercog(Uidx(ii));
    end
    return;
end

% =========single rec from here==============

DBFN = fullfile(processed_datadir, 'JL', 'Jercog', [mfilename '.dbase']);
init_dbase(DBFN, 'UniqueRecordingIndex', 'onlyifnew');
CFN = fullfile(processed_datadir, 'JL', 'Jercog', mfilename);
[J, CFN, CP] = getcache(CFN, Uidx);
if isempty(J), % compute
    db = JLdbase(Uidx);
    W = JLwaveforms(Uidx);
    WS = JLcycleStats(Uidx);
    FN1 = {'VarSpont1'    'VarOnset'    'VarDriv'  'VarTail'    'VarOffset'    'VarSpont2'};
    FN2 = {'MaxMeanIpsiCycle'  'MinMeanIpsiCycle'  'MaxVarIpsiCycle'  'MinVarIpsiCycle' ...
        'MaxMeanContraCycle'    'MinMeanContraCycle'    'MaxVarContraCycle' ...
        'MinVarContraCycle'    'P2PMeanIpsiCycle'    'P2PVarIpsiCycle'    'P2PMeanContraCycle'    ...
        'P2PVarContraCycle'};
    maxSlope_I = max(diff(W.IpsiMeanrec)/W.dt_IpsiMean);
    maxSlope_C = max(diff(W.ContraMeanrec)/W.dt_ContraMean);
    J = structJoin(db, '-Jercog', CollectInStruct(maxSlope_I, maxSlope_C));
    J = structJoin(J, '-spont', structpart(WS, FN1), '-Peax', structpart(WS, FN2));
    putcache(CFN, 1e4, CP, J);
    Add2dbase(DBFN, J);
end



