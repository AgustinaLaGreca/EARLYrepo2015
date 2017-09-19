function SP = JLspikes(Uidx, doPlot)
% JLspikes - view spikes timing and spike waveforms of JL data
%    SP = JLspikes(Uidx, doPlot)

[doPlot] = arginDefaults('doPlot', nargout<1);

%S = JLdatastruct(Uidx); % standardized recording info
D = JLdbase(Uidx); % general info on the recording & stimulus
S = JLdatastruct(Uidx); % specific info on stimulus etc
W = JLwaveforms(Uidx); % waveforms and duration params
T = JLcycleStats(Uidx); % cycle-based means & variances
A = JL_APana(Uidx,0); % EPSP-AP latencies
AD = JLanovaDetails(Uidx);
TT = JLtitle(Uidx);
dt = W.dt;

if ~isfield(T, 'VarTail'), % get new versions of T & W
    [T W] = JLviewrec(Uidx, 0); % 0 = don't plot
    '======================='
end

SPT = W.SPTraw(:, W.APinStim)-T.t_stimOffset;
if isempty(SPT), SPT = []; end; % avoid pedantic errors re sizes
lat_st = A.Lat;
if isempty(lat_st), lat_st = []; end; % avoid pedantic errors re sizes
EPSPT = SPT-lat_st;
Nspike = numel(SPT);
Tbeat = 1e3/S.Fbeat; % beat period in ms
T_I = 1e3/S.Freq1; % ipsi cycle in ms
T_C = 1e3/S.Freq2; % contra cycle in ms
SPT_I = mod(SPT, T_I)/T_I;
SPT_C = mod(SPT, T_C)/T_C;
SPT_B = mod(SPT, Tbeat)/Tbeat;
[VS_I, vs_alpha_I] = vectorstrength(SPT, S.Freq1);
[VS_C, vs_alpha_C] = vectorstrength(SPT, S.Freq2);
[VS_B, vs_alpha_B] = vectorstrength(SPT, S.Fbeat);
EPSPT_I = mod(EPSPT, T_I)/T_I;
EPSPT_C = mod(EPSPT, T_C)/T_C;
EPSPT_B = mod(EPSPT, Tbeat)/Tbeat;
[eVS_I, evs_alpha_I] = vectorstrength(EPSPT, S.Freq1);
[eVS_C, evs_alpha_C] = vectorstrength(EPSPT, S.Freq2);
[eVS_B, evs_alpha_B] = vectorstrength(EPSPT, S.Fbeat);

UniqueRecordingIndex = Uidx;
SP = CollectInStruct(UniqueRecordingIndex, ...
    '-vector_strength', VS_I, VS_C, VS_B, vs_alpha_I, vs_alpha_C, vs_alpha_B, ...
    eVS_I, eVS_C, eVS_B, evs_alpha_I, evs_alpha_C, evs_alpha_B);
GBdir = fullfile(processed_datadir,  '\JL\NNTP');
DBFN = fullfile(GBdir, 'PhaseLock.dbase');
init_dbase(DBFN, 'UniqueRecordingIndex', 'onlyifnew');
Add2dbase(DBFN, SP);

if ~doPlot, return; end
%===========snippets==========
set(figure,'units', 'normalized', 'position', [-0.0188 0.562 0.333 0.354])
LegStr = {}; H = [];
if any(W.APinSpont1), 
    h = plot(W.Tsnip, W.Snip(:, W.APinSpont1), 'k'); 
    H = [H, h(1)];
    LegStr = [LegStr [num2str(T.NspikeInSpont1) ' spont1']];
end
if any(W.APinOnset), 
    h = xplot(W.Tsnip, W.Snip(:, W.APinOnset), 'r'); 
    H = [H, h(1)];
    LegStr = [LegStr [num2str(T.NspikeInOnset) ' onset']];
end;
if any(W.APinStim), 
    h = xplot(W.Tsnip, W.Snip(:, W.APinStim), 'linewidth', 2); 
    H = [H, h(1)];
    LegStr = [LegStr [num2str(T.NspikeInStim) ' steady state']];
end;
if any(W.APinTail), 
    h = xplot(W.Tsnip, W.Snip(:, W.APinTail), 'b'); 
    H = [H, h(1)];
    LegStr = [LegStr [num2str(T.NspikeInTail) ' tail']];
end;
if any(W.APinOffset), 
    h = xplot(W.Tsnip, W.Snip(:, W.APinOffset), 'g'); 
    H = [H, h(1)];
    LegStr = [LegStr [num2str(T.NspikeInOffset) ' offset']];
end;
if any(W.APinSpont2), 
    h = xplot(W.Tsnip, W.Snip(:, W.APinSpont2), 'color', [0 0.6 0]); 
    H = [H, h(1)];
    LegStr = [LegStr [num2str(T.NspikeInSpont2) ' spont 2']];
end;
xlabel('Time (ms)', 'fontsize', 14);
ylabel('Potential (mV)', 'fontsize', 14)
legend(H, LegStr)
title(TT);

%===========cycle histograms: APs==========
set(figure,'units', 'normalized', 'position', [0.313 0.554 0.17 0.372])
if Nspike>0,
    Nhist = min(50,round(Nspike/7));
    Nhist = max(Nhist, 10);
    Ph = ((1:Nhist)-0.5)/Nhist;
    ah1=subplot(3,1,1);
    hist(SPT_I, Ph);
    set(findobj(ah1,'type', 'patch'), 'FaceColor', [0 0 0.75])
    text(0.3, 0.8, ['ipsi: R=' num2str(abs(VS_I),2)], 'units', 'normalized');
    YL1 = ylim;
    title(TT);
    ah2=subplot(3,1,2);
    hist(SPT_C, Ph);
    set(findobj(ah2,'type', 'patch'), 'FaceColor', [0.75 0 0])
    text(0.3, 0.8, ['contra: R=' num2str(abs(VS_C),2)], 'units', 'normalized');
    YL2 = ylim;
    ah3=subplot(3,1,3);
    hist(SPT_B, Ph);
    set(findobj(ah3,'type', 'patch'), 'FaceColor', [0 0.75 0])
    text(0.3, 0.8, ['bin: R=' num2str(abs(VS_B),2)], 'units', 'normalized');
    xlabel('Phase (cycle)');
    YL3 = ylim;
    H = [ah1 ah2 ah3];
    set(H, 'ylim' ,[0 1.2*max([YL1, YL2, YL3])]);
else,
    text(0.5,0.5, 'no steady-state spikes', 'horizontalalign', 'center')
end

%===========cycle histograms: EPSPs==========
set(figure,'units', 'normalized', 'position', [0.484 0.553 0.17 0.37])
if Nspike>0,
    ah1=subplot(3,1,1);
    hist(EPSPT_I, Ph);
    set(findobj(ah1,'type', 'patch'), 'FaceColor', [0 0 0.75])
    text(0.3, 0.8, ['ipsi: Re=' num2str(abs(eVS_I),2)], 'units', 'normalized');
    YL1 = ylim;
    title(TT);
    ah2=subplot(3,1,2);
    hist(EPSPT_C, Ph);
    set(findobj(ah2,'type', 'patch'), 'FaceColor', [0.75 0 0])
    text(0.3, 0.8, ['contra: Re=' num2str(abs(eVS_C),2)], 'units', 'normalized');
    YL2 = ylim;
    ah3=subplot(3,1,3);
    hist(EPSPT_B, Ph);
    set(findobj(ah3,'type', 'patch'), 'FaceColor', [0 0.75 0])
    text(0.3, 0.8, ['bin: Re=' num2str(abs(eVS_B),2)], 'units', 'normalized');
    xlabel('Phase (cycle)');
    YL3 = ylim;
    H = [ah1 ah2 ah3];
    set(H, 'ylim' ,[0 1.2*max([YL1, YL2, YL3])]);
else,
    text(0.5,0.5, 'no steady-state spikes', 'horizontalalign', 'center')
end

% EPSP-AP latencies
set(figure,'units', 'normalized', 'position', [0.652 0.518 0.362 0.344]);
[dum, iok] = denan(lat_st);
minLat = min(lat_st); maxLat = max(lat_st);
CM = colormap; Ncol = size(CM,1);
Lat2ci = @(lat)1+floor((Ncol-1)*lat/A.max_lat);
EPSP_phase_I = mod(bsxfun(@minus, AD.Spike_phase_I, lat_st/T_I),1); % Ipsi stim phase at the time of the EPSP
EPSP_phase_C = mod(bsxfun(@minus, AD.Spike_phase_C, lat_st/T_C),1); % ditto Contra
for ii=iok, 
    lat = lat_st(ii);
    col = CM(Lat2ci(lat),:);
    xplot(EPSP_phase_C(ii), EPSP_phase_I(ii), '*', 'color', col); 
end
h = colorbar;
LatTick = 0:0.05:0.5; 
LatTickLabel = Words2cell(num2str(round(1000*LatTick)));
set(h,'ytick', Lat2ci(LatTick), 'yticklabel', LatTickLabel);
ylabel(h, 'EPSP-AP latency (\mus)');
ylabel('Ipsi phase (cycle)');
xlabel('Contra phase (cycle)');













