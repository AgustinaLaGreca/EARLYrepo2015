function [C, SV] = JL_compareSpikeCounts(Nbin);
% JL_compareSpikeCounts - compare spike counts btw monaural & binaural stim
%   [C, SV] = JL_compareSpikeCounts;

Nbin = arginDefaults('Nbin',10);

SV = retrieve_dbase('D:\processed_data\JL\JLsplitvar\JLsplitvar.dbase');
DBFN = ['D:\processed_data\JL\JLcompareSpikeCounts\JLcompareSpikeCounts_' num2str(Nbin) '.dbase'];
init_dbase(DBFN, 'UniqueTripletIndex', 'onlyifnew')
for ii=1:numel(SV),
    if rem(ii,10)==0, disp(ii); end
    sv = SV(ii);
    unpack(sv);
    CS_I = JLcycleStats(sv.Uidx_I);
    CS_C = JLcycleStats(sv.Uidx_C);
    CS_B = JLcycleStats(sv.Uidx_B);
    %I
    NspikeInSpont1_I = CS_I.NspikeInSpont1;
    NspikeInStim_I = CS_I.NspikeInStim;
    DrivenRate_I = CS_I.DrivenRate;
    SpontRate1_I = CS_I.SpontRate1;
    %C
    NspikeInSpont1_C = CS_C.NspikeInSpont1;
    NspikeInStim_C = CS_C.NspikeInStim;
    DrivenRate_C = CS_C.DrivenRate;
    SpontRate1_C = CS_C.SpontRate1;
    %B
    NspikeInSpont1_B = CS_B.NspikeInSpont1;
    NspikeInStim_B = CS_B.NspikeInStim;
    DrivenRate_B = CS_B.DrivenRate;
    SpontRate1_B = CS_B.SpontRate1;
    bh = local_binhist(sv.Uidx_B, Nbin);
    c = CollectInStruct(...
        UniqueTripletIndex, Uidx_I, Uidx_C, Uidx_B, ...
        '-ID', iexp, icell, icell_run, iseries, iseries_run, freq, SPL, ...
        '-SpontSpikes', ...
        NspikeInSpont1_I, NspikeInStim_I, DrivenRate_I, SpontRate1_I, ...
        NspikeInSpont1_C, NspikeInStim_C, DrivenRate_C, SpontRate1_C, ...
        NspikeInSpont1_B, NspikeInStim_B, DrivenRate_B, SpontRate1_B);
    % append spikestats 
    spst = JLspikeStats(sv.Uidx_B);
    spst = rmfield(spst,'UniqueRecordingIndex');
    C(ii)= structJoin(c, bh ,spst);
end
Add2dbase(DBFN, C);


function S = local_binhist(Uidx, Nbin);
BinEdges = linspace(0,1,Nbin+1);
%D = JLdbase(Uidx); % general info on the recording & stimulus
S = JLdatastruct(Uidx); % specific info on stimulus etc
W = JLwaveforms(Uidx); % waveforms and duration params
T = JLcycleStats(Uidx); % cycle-based means & variances
dt = W.dt;
%W T S
SPT = W.SPTraw(:, W.APinStim)-T.t_stimOffset;
if isempty(SPT), SPT = []; end; % avoid pedantic errors re sizes
Nspike = numel(SPT);
Tbeat = 1e3/S.Fbeat; % beat period in ms
SPT_B = mod(SPT, Tbeat)/Tbeat;
Nsp = histc(SPT_B, BinEdges);
if isempty(Nsp),
    Nsp = zeros(Nbin+1,1);
end
Nsp = Nsp(1:Nbin); % discard garbage bin
Nsp_best = max(Nsp);
Nsp_zero = mean(Nsp([1 end]));
AnaDur = T.AnaDur;
Rate_B_best = 1e3*Nsp_best/AnaDur*Nbin; % spike rate @ best ITD, accounting for reduced duration due to binning
Rate_B_zero = 1e3*Nsp_zero/AnaDur*Nbin; % spike rate @ best ITD, accounting for reduced duration due to binning
S = CollectInStruct('-Bin_ITD', Nsp, AnaDur, Nsp_best, Nsp_zero, Rate_B_best, Rate_B_zero);
