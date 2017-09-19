function [S, mWvI, mWvC] = JLpeakPhase(Uidx, cdelay, doPlot);
% JLpeakPhase - stimulus phase at peak locations
%    JLpeakPhase(Uidx, cdelay) returns a struct having fields
%    UniqueRecordingIndex: equal to Uidx
%            cdelay: equal to input cdelay, i.e. delay compensation in ms.
%            T_ipsi: ipsi stimulus period (ms)
%         Freq_ipsi: ipsi frequency (Hz)
%          T_contra: contra stimulus period (ms)
%        Freq_contra: contra frequency (Hz)
%      Ph_peak_ipsi: phases (cycles) of peaks in ipsi-averaged waveform
%    Ph_peak_contra: phases (cycles) of peaks in ipsi-averaged waveform
%   All phases are compensated for a delay of cdelay ms. 
%
%    JLpeakPhase(-Uidx) returns struct array whose elements correspond to
%    the individual recordings of the measurement series of which Uidx is a
%    member. It also generates contour plots, etc.
%
%    JLpeakPhase(..., cdelay, doPlot) enables or disables plotting
%    depending on logical doPlot. Default is plot only if Uidx is negative
%    and no output args are requested.

[cdelay, doPlot] = arginDefaults('cdelay/doPlot', 0, nargout<1 && (numel(Uidx)>1 || Uidx<0));
if numel(Uidx)>1,
    for ii=1:numel(Uidx),
        [S(ii), mWvI(:,ii), mWvC(:,ii)] = JLpeakPhase(Uidx(ii), cdelay);
    end
    if doPlot,
        figh = local_plot(S, mWvI, mWvC);
        setGUIdata(figh, 'Uidx', Uidx);
        setGUIdata(figh, 'cdelay', cdelay);
    end
    return
elseif Uidx<0, % find series to which abs(Uidx) belongs
    db = JLdbase(abs(Uidx));
    db = JLdbase('iseries_run', db.iseries_run);
    [S, mWvI, mWvC] = JLpeakPhase([db.UniqueRecordingIndex], cdelay, doPlot);
    return;
end

W = JLwaveforms(Uidx);
D = JLdatastruct(Uidx);
Freq_ipsi = D.Freq1; % Hz
Freq_contra = D.Freq1; % Hz
T_ipsi = 1e3/Freq_ipsi; % ms
T_contra = 1e3/Freq_contra; % ms

[Ph_peak_ipsi, mWvI]= local_peak(W, T_ipsi, 'Ipsi', cdelay);
[Ph_peak_contra, mWvC] = local_peak(W, T_contra, 'Contra', cdelay);
ShiftAddPeak = local_SAP(mWvI, mWvC);

UniqueRecordingIndex = Uidx;
S = CollectInStruct(UniqueRecordingIndex, '-WvParams', cdelay, T_ipsi, Freq_ipsi, ...
    T_contra, Freq_contra, '-peakphases', Ph_peak_ipsi, Ph_peak_contra, ShiftAddPeak);

%============================================================

function [Ph_peak, mWv] = local_peak(W, T, Ch, cdelay);
if isnan(cdelay), cdelay=0; end
CS = JLcycleStats(W.UniqueRecordingIndex);
mWv = W.([Ch 'Meanrec']); % subtracting  CS.MeanSpont1 gave unwanted effects
mWv = mWv - mean(mWv);
% resample to get uniform sizes
NsamOld = numel(mWv);
NsamNew = 500;
mWv= interp1(linspace(0,1,NsamOld)', mWv, linspace(0,1,NsamNew)');
dph = 1/NsamNew;
% apply time shift realizing cdelay
dt = W.(['dt_' Ch 'Mean'])*NsamOld/NsamNew;
NsamShift = round(-cdelay/dt);
mWv = circshift(mWv, NsamShift);
Ph_peak = peakfinder(dph, repmat(mWv,[2 1]), 0.3/T); % take 2 periods to eliminate boundary maxima
Ph_peak = Ph_peak(Ph_peak>0.5 & Ph_peak<1.5);
Ph_peak = sort(mod(Ph_peak,1));

function WB = local_SAP(I, C);
Im = cconvmtx(I(:)); % columns contain all circ-shifted versions of I
C = repmat(C(:), [1 size(Im,2)]);
WB = max(Im+C).';


function figh = local_plot(S, mWvI, mWvC);
set(figure,'units', 'normalized', 'position', [0.0438 0.279 0.891 0.621], 'toolbar', 'figure');
figh = gcf;
XY = [850 600];
CLR = get(figh,'color');
UIp = struct('fontsize', 12);
uicontrol(figh, 'position', [XY 80 20], 'style', 'text', 'string', 'cdelay:', 'backgroundcolor', CLR, UIp);
he=uicontrol(figh, 'position', [XY+[70 0] 50 20], 'style', 'edit', 'string', num2str(S(1).cdelay), 'backgroundcolor', [1 1 1], 'horizontalalign', 'left', UIp);
uicontrol(figh, 'position', [XY+[110 0] 40 20], 'style', 'text', 'string', 'ms', 'backgroundcolor', CLR, UIp);
uicontrol(figh, 'position', [XY+[150 0] 60 20], 'string', 'submit', 'backgroundcolor', CLR, UIp, 'callback', @(Src,Ev)local_subm(Src,he));
subplot(2,2,1);
for ii=1:numel(S), 
    fr = S(ii).Freq_ipsi; 
    ph = S(ii).Ph_peak_ipsi; 
    ph = mod(ph+0.5,1)-0.5;
    xplot(fr+0*ph, ph, 'b*'); 
    ph = S(ii).Ph_peak_contra; 
    ph = mod(ph+0.5,1)-0.5;
    xplot(fr+0*ph, ph, 'rs'); 
end
ylim([-0.5 0.5]);
ylabel('Monaural phase (cycles)', 'fontsize',14);
DB = JLdbase([S.UniqueRecordingIndex]);
db = DB(1);
TT = sprintf('Exp %d cell %d (icell_run=%d)  %s %d dB  series %d', db.iexp, db.icell, db.icell_run, db.chan, db.SPL, db.iseries);
title([blanks(70) TT], 'interpreter', 'none', 'fontsize', 16);
%=
subplot(2,2,3);
Nph = size(mWvI,1); % # phase values
Nh = round(Nph/2);
dph = 1/Nph;
PH = dph*(0:Nph-1).'-0.5;
[Wmin, Wmax] = minmax([mWvI(:); mWvC(:)]);
contourf([S.Freq_ipsi], PH, circshift(mWvI ,Nh));
caxis([Wmin, Wmax]);
xlabel('Frequency (Hz)', 'fontsize',14);
ylabel('Monaural phase (cycles)', 'fontsize',14);
title('IPSI', 'fontsize',12)
%=
subplot(2,2,4);
contourf([S.Freq_ipsi], PH, circshift(mWvC ,Nh));
caxis([Wmin, Wmax]);
title('CONTRA', 'fontsize',12)
xlabel('Frequency (Hz)', 'fontsize',14);
ylabel('Monaural phase (cycles)', 'fontsize',14);
%=
subplot(2,2,2);
SAP = [S.ShiftAddPeak];
Freq = [S.Freq_ipsi];
contourf(Freq, PH, circshift(SAP ,Nh));
ylabel('Binaural phase (cycles)', 'fontsize',14);
caxis(max(SAP(:))+[2*(Wmin-Wmax) 0]);
for ii=1:numel(S),
    Uidx = S(ii).UniqueRecordingIndex;
    D = JLdatastruct(Uidx);
    SPT = JLspiketimes(Uidx); % in ms since stim onset
    Tbeat = 1e3/D.Fbeat; % ms
    SPH = mod(SPT/Tbeat+0.5,1)-0.5;
    SetRandState(Uidx);
    dfreq = 10*randn(size(SPH));
    xplot(dfreq+Freq(ii)+5, SPH, '.k', 'markersize',7);
    xplot(dfreq+Freq(ii), SPH, '.w', 'markersize',7);
end


%====================================
function local_subm(Src,he)
figh = parentfigh(Src);
Uidx = getGUIdata(figh, 'Uidx');
cdelay = str2num(get(he, 'string'));
if ~isempty(cdelay),
    close(figh);
    JLpeakPhase(Uidx, cdelay);
end







