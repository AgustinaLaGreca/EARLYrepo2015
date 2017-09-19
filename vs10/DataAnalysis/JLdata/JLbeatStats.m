function ST = JLbeatStats(S, iCell, doCompute);
% JLbeatStats - evaluate stats of JLbeat analysis
%   ST = JLbeatStats(S), where S is output of JLbeat, returns struct
%   array ST containing elementary stats of S.
%
%   ST = JLbeatStats(ExpID, iCell) calls JLreadBeat to retrieve the processed
%   data of experiment ExpID, cell # iCell, and concatenates the available series
%   of beat data.
%
%   JLbeatStats(ExpID) or JLbeatStats(ExpID, 0) concatenates all stats available for 
%   experiment ExpID. Data retrieval uses JLreadBeats('list').
% 
%   JLbeatStats('all') combines data from all experiments.
%   
%   See also JLbeat, JLbeatPlot, JLreadBeats.

%     ExpID
%     RecID
%     seriesID
%     iGerbil
%     icond
%     icell
%     a___________________
%     UniqueCellIndex
%     UniqueSeriesIndex
%     UniqueRecordingIndex
%     ABFname
%     ABFdir
%     JLcomment
%     MHcomment
%     b___________________
%     Nxmean             # fine structure cycle averages used to evaluate beat waveform 
%     APthrSlope         max downward slope criterion [V/s] for APs
%     APwindow_start     start of AP window re AP peak [ms]
%     APwindow_end       end of AP window re AP peak   [ms]
%     CutoutThrSlope     threshold [V/s] for truncating peak w/o necessarily counting the peak as an AP
%     c___________________
%     StimType      'BFS' for beats
%     Freq1         ipsi freq [Hz]
%     Freq2         contra freq [Hz]
%     Fbeat         beat freq [Hz]
%     burstDur      duration [ms] of tone
%     SPL1          intensity [dB SPL] of ipsi tone
%     SPL2          intensity [dB SPL] of contra tone
%     d___________________
%     Fsam        sample rate [Hz] of recording
%     dt          sample period [ms] of recording
%     Nbeats      # beats in recording
%     AnaDur      length of analysis window [ms]
%     Ttrans      duration [ms] of transient response 
%     MaxOrd      maximum order of distortions considered in spectral analysis
%     TTT         char string identifying the stimulus condition
%     Period1     period [ms] of ipsi freq
%     Period2     period [ms] of contra freq
%     e___________________
%     f___________________
%     df          spacing [Hz] of spectrum
%     Yoffset     vertical offset for waveform plots
%     AC          shuffled autocorrelation of waveform across beat cycles
%     Tshift      time shift (ms) used to move beat max to middle of plot
%     g___________________
%     h___________________
%     varName
%     st___Partners_______
%     Ipsi_Partners      recording IDs of ipsi measurements using same stim
%     Contra_Partners    recording IDs of contra measurements using same stim
%     Bin_Partners       recording IDs of binaural measurements using same stim
%     st___Wvstats________
%     Vspont_STD      STD [mV] during spontaneous activity (first 500 ms of recording)
%     Vstim_STD       STD [mV] during the analysis window
%     Vbin_STD        STD [mV] of the beat-cycle average waveform
%     Vbin_P2P        Peak-Peak [mV] of the beat-cycle average (BCA) waveform
%     Vipsi_STD       STD [mV] of the ipsi contribution to BCA waveform
%     Vipsi_P2P       Peak-Peak [mV] of the ipsi contribution to BCA waveform
%     Vipsi_Npeak     # peaks in the ipsi contribution to BCA waveform
%     Vcontra_STD     STD [mV] of the contra contribution to BCA waveform
%     Vcontra_P2P     STD [mV] of the contra contribution to BCA waveform
%     Vcontra_Npeak   # peaks in the contra contribution to BCA waveform
%     st___WvPred_________
%     VarAccLin       % variance of BCA accounted for by adding ipsi&contra contributions
%     CorrLin         Pearson correlation btw BCA and linear pred
%     VarAccInteract  % var. of BCA acc. for by nonlinear interactin of ipsi & contra contributions
%     CorrInteract    Pearson correlation btw BCA and nonlin interaction pred
%     st___Spikes_________
%     Nspike          # spikes during analysis window
%     SpikeRate       rate (spikes/s) during analysis window
%     st___STA____________
%     STAbin_Npk      # peaks in spike-triggered average (STA) of BCA
%     STAbin_Tpk1     time of largest STA peak re AP peak
%     STAbin_P2P      Peak-Peak [mV] of STA
%     STAipsi_Npk     # peaks in STA of ipsi-cycle-averaged waveform
%     STAipsi_Tpk1    # time of largest peak in STA of ipsi-cycle-averaged waveform
%     STAipsi_P2P     Peak-PEak [mV] of STA of of ipsi-cycle-averaged waveform
%     STAcontra_Npk   # peaks in STA of contra-cycle-averaged waveform
%     STAcontra_Tpk1  # time of largest peak in STA of contra-cycle-averaged waveform
%     STAcontra_P2P   Peak-PEak [mV] of STA of of contra-cycle-averaged waveform
%     st___ph_lock________
%     VSbin           beat vector strength 
%     Raleighbin      Rayhleigh stat of binaural vector strength
%     PHbin           phase of binaural vector strength
%     VSipsi          ipsi vector strength
%     Raleighipsi     Rayhleigh stat of ipsi vector strength
%     PHipsi          phase of ipsi vector strength
%     VScontra        contra vector strength
%     Raleighcontra   Rayhleigh stat of contra vector strength
%     PHcontra        phase of contra vector strength
%     st___spont__________
%     SpontAPcount    # APs in 500-ms spontaneous part
%     SpontDTsmooth   smoothing window [ms] used for spont analyis
%     SpontPeakWidth  peak width used for analysis of spont peaks
%     SpontAcceptCrit sharpness criterium used for detecting spon peaks
%     NspontPeak      # peaks in spont trace
%     SpontEventRate  peak rate [peaks/s] in spont trace

if nargin<1, S = 'all'; end
if nargin<2, iCell=0; end
if nargin<3, doCompute=0; end % do not (re-)compute JLbeat stuff

if isequal('all', S), % 4th syntax of help text: ST = JLbeatStats('all'). Recursive handling.
    L = JLreadBeats('list');
    Elist = fieldnames(L);
    ST = [];
    for ii=1:numel(Elist),
        ExpID = Elist{ii};
        st = JLbeatStats(ExpID, 0, doCompute);
        ST = [ST, st];
    end
    return;
elseif ischar(S), % 2nd or 3rd syntax of help text: ST = JLbeatStats(ExpID, iCell). Recursive handling.
    if isequal(0,iCell), % get list of cells
        L = JLreadBeats('list');
        iCell = L.(upper(S));
    end
    ST = [];
    for ic=iCell(:).',
        ic
        JLreadBeats(S, ic, doCompute); % this produces local variables named Jb*
        dseries = who('Jb*');
        for ii=1:numel(dseries), % recursive call to process a single series
            Jb = eval([dseries{ii} ';']);
            st = JLbeatStats(Jb);
            ST = [ST, st];
        end
    end
    return;
end

% copy S into ST while removing waveforms and other array-valued fields
rflds = Words2cell('Stim1/Stim2/R0/Y1/Y2/Y3/Yc/Sb/Fdpb/SPTraw/SPT/Tsnip/Snip/APslope','/');
ST = rmfield(S, rflds);

% todo: max derivatives of I&C wv, peakiness. Think about Plotting
%    baseline RMS
%    absolute recording time (chronology)
%    systematic catalog of recrdings per cell


% compute stats
for ii=1:numel(S),
    s = S(ii);
    ST(ii).(local_sepname('Partners')) = '________________';
    try,
        PRT = JLpartners(s);
        ST(ii).Ipsi_Partners = PRT.Ipsi;
        ST(ii).Contra_Partners = PRT.Contra;
        ST(ii).Bin_Partners = PRT.Bin;
        ST(ii).(local_sepname('Wvstats')) = '________________';
    catch
        disp('identification of partners requires 2nd run of JLbeatStats');
        ST(ii).Ipsi_Partners = [];
        ST(ii).Contra_Partners = [];
        ST(ii).Bin_Partners = [];
        ST(ii).(local_sepname('Wvstats')) = '________________';
    end
    % ==R0: beat-averaged recording
    [ST(ii).Vspont_STD, ST(ii).Vstim_STD] = local_raw_STD(s);
    ST(ii).Vbin_STD = std(s.R0);
    ST(ii).Vbin_P2P = max(s.R0)-min(s.R0);
    % ==Yl: ipsi-averaged recording
    ST(ii).Vipsi_STD = std(s.Y1);
    ST(ii).Vipsi_P2P = max(s.Y1)-min(s.Y1);
    ST(ii).Vipsi_Npeak = local_Npeak(s,1);
    % ==Y2: contra-averaged recording
    ST(ii).Vcontra_STD = std(s.Y2);
    ST(ii).Vcontra_P2P = max(s.Y2)-min(s.Y2);
    ST(ii).Vcontra_Npeak = local_Npeak(s,2);
    % ==Y3 and Yc, the predicted binaural responses
    ST(ii).(local_sepname('WvPred')) = '________________';
    ST(ii).VarAccLin = 100*(1-var(s.R0-s.Yc)/var(s.R0));
    ST(ii).CorrLin = corr(s.R0, s.Yc);
    ST(ii).VarAccInteract = 100*(1-var(s.R0-s.Y3)/var(s.R0));
    ST(ii).CorrInteract = corr(s.R0, s.Y3);
    % generic spike stats
    ST(ii).(local_sepname('Spikes')) = '________________';
    ST(ii).Nspike = numel(s.SPT);
    ST(ii).SpikeRate = 1e3*numel(s.SPT)/s.AnaDur; % sp/s
    % spike-triggered averages
    ST(ii).(local_sepname('STA')) = '________________';
    TimeWin =  [-2 1];
    [ST(ii).STAbin_Npk,  ST(ii).STAbin_Tpk1, ST(ii).STAbin_P2P] = local_STA(s.dt, s.R0, s.SPT, TimeWin);
    [ST(ii).STAipsi_Npk, ST(ii).STAipsi_Tpk1, ST(ii).STAipsi_P2P] = local_STA(s.dt, s.Y1, s.SPT, TimeWin);
    [ST(ii).STAcontra_Npk, ST(ii).STAcontra_Tpk1, ST(ii).STAcontra_P2P] = local_STA(s.dt, s.Y2, s.SPT, TimeWin);
    % phase locking
    ST(ii).(local_sepname('ph_lock')) = '________________';
    [R1, alpha1] = vectorstrength(s.SPT,s.Freq1);
    [R2, alpha2] = vectorstrength(s.SPT,s.Freq2);
    [Rb ,alphab] = vectorstrength(s.SPT,s.Fbeat);
    ST(ii).VSbin = abs(Rb);
    ST(ii).Raleighbin = alphab;
    ST(ii).PHbin = cangle(Rb);
    ST(ii).VSipsi = abs(R1);
    ST(ii).Raleighipsi = alpha1;
    ST(ii).PHipsi = cangle(R1);
    ST(ii).VScontra = abs(R2);
    ST(ii).Raleighcontra = alpha2;
    ST(ii).PHcontra = cangle(R2);
    % spontaneous stuff
    ST(ii).(local_sepname('spont')) = '________________';
    spst = JLspontPeaks(s,0.3, []);
    ST(ii).SpontAPcount = spst.NspontAP;
    ST(ii).SpontDTsmooth = spst.DTsmooth;
    ST(ii).SpontPeakWidth = spst.PeakWidth;
    ST(ii).SpontAcceptCrit = spst.AcceptCrit;
    ST(ii).NspontPeak = spst.Npeak;
    ST(ii).SpontEventRate = spst.SpontEventRate;
    ST(ii).SpontMedianDVmin = median(spst.DVmin);
    ST(ii).SpontMedianDVmax = median(spst.DVmax);
end



%=============================
function n = local_sepname(c); 
% add separator
n = ['st___' c repmat('_', [1 15-numel(c)])];

function [Npk, Tpk1, P2P] = local_STA(dt, Wv, SPT, TimeWin);
STA = SpikeTrigAv(dt, Wv, SPT, TimeWin);
Tw = timeaxis(STA, dt, TimeWin(1));
[Tpk, Ypk, isort]=localmax(Tw,STA ,0.15);
Npk = numel(Tpk);
if Npk>0,
    [Tpk, Ypk] = deal(Tpk(isort), Ypk(isort));
    Tpk1 = Tpk(1);
else,
    Tpk1 = nan;
end
P2P = max(STA)-min(STA);

function Npeak = local_Npeak(s,iear);
if isequal(1,iear),
    y = smoothen(s.Y1,0.3, s.dt);
    T = 1e3/s.Freq1; % stim period in ms
else,
    y = smoothen(s.Y2,0.3, s.dt);
    T = 1e3/s.Freq2; % stim period in ms
end
Nsam = numel(y);
Ncycle = Nsam*s.dt/T; % # stim cycles in y
Tpeak = localmax(s.dt,y);
%f3 ;dplot(s.dt, y); fenceplot(Tpeak, ylim);
Npeak = round(numel(Tpeak)/Ncycle); % # peaks per stim cycle

function  [STD_spont, STD_stim] = local_raw_STD(s);
StimOnset = 500; % ms start of beat stim
HPfreq = 20; % Hz highpass freq
D = readTKABF(s);
dt = D.dt_ms;
R0 = D.AD(1).samples; clear D;
R0 = APtruncate2(R0, s.CutoutThrSlope, dt, [s.APwindow_start s.APwindow_end]);
[B,A] = butter(5, 2e-3*HPfreq*dt ,'high');
R = filter(B,A,R0);
%dplot(dt, R); xdplot(dt, R0, 'r'); pause; delete(gcf);
NsamSpont = 1+round(StimOnset/dt);
STD_spont = std(R(1:NsamSpont));
AnaOnset = StimOnset  + s.Ttrans;
AnaOffset = AnaOnset + s.AnaDur;
isam0 = 1+round(AnaOnset/dt);
isam1 = round(AnaOffset/dt);
STD_stim = std(R(isam0:isam1));




