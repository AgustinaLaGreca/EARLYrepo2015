function S = JLbestITD(Uidx, Add2DB);
% JLbestITD - best ITD in us
%    S = JLbestITD(Uidx)
%    S is struct (array) holding Uidx, bestITD, freq, SPL.

Add2DB = arginDefaults('Add2DB', 0);

if numel(Uidx)>1,
    for ii=1:numel(Uidx),
        S(ii) = JLbestITD(Uidx(ii));
    end
    return;
end

%======= single Uidx from here=============
Sidx = floor(Uidx/100);

CFN = fullfile(processed_datadir, 'JL\JLbestITD\bestITD');
[sITD, CFN, CP] = getcache(CFN, Sidx);
if isempty(sITD),
    clear sITD
    D = JLdatastruct(Uidx);
    DB = JLdbase('iseries_run', D.iseries_run); % all recs of this series
    for ii=1:numel(DB),
        sITD(ii) = local_doit(DB(ii).UniqueRecordingIndex);
    end
    sITD = sortAccord(sITD, [sITD.freq]);
    putcache(CFN, 5e3, CP ,sITD);
end
S = sITD([sITD.UniqueRecordingIndex]==Uidx);
if Add2DB,
    DBFN = fullfile(processed_datadir, 'JL\JLbestITD\bestITD.dbase');
    init_dbase(DBFN, 'UniqueRecordingIndex', 'onlyifnew');
    Add2dbase(DBFN, S);
end

function S = local_doit(Uidx);
Nbin = 12;
Edge = (0:Nbin)/Nbin-0.5;
BC = (Edge(2:end)+Edge(1:end-1))/2;
D = JLdatastruct(Uidx);
ST = JLspikeStats(Uidx);
SPT = JLspiketimes(Uidx);
S_ph = mod(SPT/250+0.5,1)-0.5;
Nsp = histc(S_ph, Edge);
if isempty(Nsp), Nsp = 0*Edge; end
Nsp(end) = [];
freq = D.Freq1;
if ~isequal('B', D.chan) || all(Nsp==0) || (ST.vs_alpha_B>0.001), % not binaural, no spikes or not binaurally phase locking
    bestITD = nan;
else,
    [dum, imax] = max(Nsp);
    bestITD = 1e6*BC(imax)/(D.Freq1+2);
end
UniqueRecordingIndex = Uidx;
SPL = D.SPL;
VS_B = ST.VS_B;
vs_alpha_B = ST.vs_alpha_B;
S = CollectInStruct(UniqueRecordingIndex, freq, SPL, bestITD, VS_B, vs_alpha_B);









