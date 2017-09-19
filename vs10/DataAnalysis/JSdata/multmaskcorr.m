function S = multmaskcorr(DS);
% multmaskcorr - analyze group of MASK datasets together
MaxLag = 10;


Nd = numel(DS);
for ii=1:Nd,
    Stim(ii) = DS(ii).Stim;
    iRec(ii) = DS(ii).ID.iDataset;
    Uidx(ii) = datenum(DS(ii).ID.created);
end

[S, CFN, CP] = getcache(fullfile(processed_datadir, 'JS\multmask', mfilename), Uidx);
if isempty(S), % compute it
    S = local_compute(DS, Stim, iRec, Uidx, MaxLag);
    putcache(CFN, 200, CP, S);
end

figure;
for iw=1:S.Nw,
    xplot(S.Tau, S.NSac{iw}, lico(iw));
end

%=================================================================
function S = local_compute(DS, Stim, iRec, Uidx, MaxLag);
Nd = numel(DS);
if numel(unique([Stim.MnoiseSPL]))>1,
    error('Non matching Noise SPLs');
end
if numel(unique([Stim.ToneFreq]))>1,
    error('Non matching tone frequencies');
end
if numel(unique([Stim.StartWarp]))>1 || numel(unique([Stim.StepWarp]))>1 || numel(unique([Stim.EndWarp]))>1,
    error('Non matching warp params.');
end
[DS, Stim] = sortAccord(DS, Stim, [Stim.ToneSPL]);
ToneSPL = [Stim.ToneSPL];
Warp = Stim(1).Warp;
Nw = numel(Warp);
Nrep = Stim(1).Nrep;
StimDur = Stim(1).MnoiseBurstDur;
ToneOnset = Stim(1).ToneOnset;
NoiseDur = ToneOnset;
ToneDur = StimDur - ToneOnset;
ToneFreq = Stim(1).ToneFreq;
NoiseSPL = Stim(1).MnoiseSPL;

if ToneFreq<=1000,
    BinWidth = 1e3/ToneFreq/12; % 1/12 of a cycle
else,
    BinWidth = 0.1;
end
% get spike times separately for noise & tone portions. Note that Anwin unwarps.
SPTn = cell(Nd, Nw, Nrep);
SPTt = cell(Nd, Nw, Nrep);
for id=1:Nd,
    for iw=1:Nw,
        SPTn(id,:,:) = AnWin(DS(id),[0 ToneOnset]);
        SPTt(id,:,:) = AnWin(DS(id),[ToneOnset StimDur]);
    end
end
% noise correlogram: pool datasets
for iw=1:Nw,
    spt = SPTn(:,iw,:);
    [NSacVar{iw} Tau NSac{iw}] = SPTcorrvar(spt(:), 'nodiag', MaxLag, BinWidth, NoiseDur, 'DriesNorm');
    for id=1:Nd,
        [TSacVar{id,iw} Tau TSac{id,iw}] = SPTcorrvar(SPTt(id,iw,:), 'nodiag', MaxLag, BinWidth, NoiseDur, 'DriesNorm');
    end
end
ExpName = expname(DS(1));
S = CollectInStruct(ExpName, iRec, Uidx, ToneSPL, Warp, StimDur, ...
    Nd, Nw, Nrep, ToneOnset, NoiseDur, ToneDur, ToneFreq, NoiseSPL, ...
    '-Corr', MaxLag, BinWidth, Tau, NSacVar, NSac, TSacVar, TSac);



















