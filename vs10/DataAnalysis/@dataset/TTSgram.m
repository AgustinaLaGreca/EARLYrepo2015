function S = TTSgram(D, irecChan, scaleFac, iselCond, doPlot);
% dataset/TTSgram - periodogram of TTS data
%   TTSgram(D, irecChan, scaleFac, iselCond, doPlot);


[scaleFac, iselCond, doPlot] = arginDefaults('scaleFac,iselCond,doPlot',1,0,true); % default: no scaling of suppressor; all conditions; do plot

if isequal(0,scaleFac),
    doScale = 0;
    scaleFac = 1;
else,
    doScale = 1;
end
St = D.Stim;
Fcutoff = min(St.Fprobe)-10*St.SupFreq; % filter freq to separate bias tone & probe-related cmps

% reference dataset
[RefNosup, RefLargeSup] = deal(dataset(), D(1)); % void and D istself, resp
if numel(D)==2, % second element is no-supp ref
    RefNosup = D(2);
elseif numel(D)==3, % second element is no-supp ref; 3rd is high-SPL supp ref
    RefNosup = D(2);
    RefLargeSup = D(3);
end
D = D(1);

GP = GenericStimparams(D);
[Scal, Yunit] = conversionfactor(anachan(D, irecChan));
if isequal('mm/s', Yunit),
    Scal = 1e6*Scal;
    DisplacementUnit = 'nm';
else,
    DisplacementUnit = 'AU';
end

S = struct('Expname', expname(D), 'irec', irec(D), ...
    'irecRefNoSup', irec(RefNosup), 'irecRefLargeSup', irec(RefLargeSup), ...
    'SupSPL', St.SupSPL, 'Fsup', St.Fsup_exact);
% high_intensity suppressor ref
[lrecRef, S.SPLRefLargeSup] = local_SupRef(D, RefLargeSup, irecChan, Fcutoff, Scal);
% title string
Tstring = ['  ', expname(D) '    irec=' num2str(irec(D)) ...
    '  supp: ' num2str(St.SupSPL) ' dB, ', num2str(round(St.Fsup_exact)) ' Hz' ...
    '  Large-SPL ref: ' num2str(S.SPLRefLargeSup) 'dB    supp scaling: ' num2str(scaleFac) ''];

if doPlot,
    fh1 = figure;
    set(fh1,'units', 'normalized', 'position', [0.214 0.09 0.741 0.792]);
    [hax, Lh, Bh] = plotpanes(GP.Ncond);
    fh2 = figure;
    set(fh2,'units', 'normalized', 'position', [0.204 0.11 0.741 0.792]);
    [hax2, Lh2, Bh2] = plotpanes(GP.Ncond);
end
MAL = 0;
Ncond = GP.Ncond;
iselCond = replaceMatch(iselCond(:).', 0, 1:Ncond);
[Lrec, Hrec, cycleDur, dt] = local_Dsplit(D, irecChan, Fcutoff);
for icond = iselCond,
    %[rec, dt] = anamean(D, irecChan, icond);
    %[lrec, hrec, cycleDur] = local_split(rec, dt, icond, Fcutoff, GP, St);
    lrec = Lrec(:,icond);
    hrec = Hrec(:,icond);
    if ~isempty(lrecRef),
        lrec = lrecRef;
    end
    MAL = max(MAL, max(abs(lrec)));
    % return args
    S.dt = dt;
    S.cycleDur = cycleDur;
    S.Fprobe(1,icond) = St.Fprobe(icond);
    S.lrec(:,icond) = lrec(:);
    S.hrec(:,icond) = hrec(:);
    delay = local_delay(real(lrec),abs(hrec),dt);
    S.delay(icond) = delay;
    if doPlot,    
        %==plot waveforms
        figure(fh1);
        axes(hax(icond));
        dplot(dt,real(hrec),'b');
        xdplot(dt,abs(hrec),'k', 'linewidth',2);
        xdplot(dt,-abs(hrec),'k', 'linewidth',2);
        xdplot(dt,scaleFac*real(lrec),'r');
        xlim([0,cycleDur]);
        text(0.1, 0.9, CondLabel(D,icond), 'units', 'normalized');
        %==plot envelop2 lissajous
        figure(fh2);
        axes(hax2(icond));
        plot(scaleFac*real(lrec),abs(hrec));
        isamqt = round(numel(hrec)/4);
        xplot(scaleFac*real(lrec(1)),abs(hrec(1)), '*')
        xplot(scaleFac*real(lrec(isamqt)),abs(hrec(isamqt)), 'o')
        text(0.1, 0.9, CondLabel(D,icond), 'units', 'normalized');
        text(0.9, 0.9, ['\tau=' num2str(delay) ' \mus'], 'units', 'normalized', 'horizontalAlign', 'right');
        if ismember(hax2(icond), Lh2),
            ylabel('Probe Envelope (nm)');
        end
        if ismember(hax2(icond), Bh2),
            xlabel('Suppressor displacement (nm)');
        end
    end
end
S.Tstring = Tstring;
S.iselCond = iselCond;

if doPlot,
    figure(fh1);
    local_Yscale(doScale, scaleFac, hax, GP.Ncond, MAL);
    local_labels(hax, Bh, Lh, D, St, DisplacementUnit, scaleFac);
    axes(hax(1));

    title(Tstring, 'fontsize', 16, 'fontweight', 'bold', 'HorizontalAlignment', 'left');
    if exist('hax2', 'var'),
        axes(hax2(1));
        title(Tstring, 'fontsize', 16, 'fontweight', 'bold', 'HorizontalAlignment', 'left');
    end
end
if ~isvoid(RefNosup), % indicate normalization
    [Lrec, Hrec, cycleDur, dt] = local_Dsplit(RefNosup, irecChan, Fcutoff);
    for icond = iselCond,
        lrec = Lrec(:,icond);
        hrec = Hrec(:,icond);
        %[lrec, hrec, cycleDur] = local_split(rec, dt, icond, Fcutoff, GP, St);
        RefAmp = rms(abs(hrec));
        % return args
        S.hrecNoSup(:,icond) = hrec(:);
        S.rmsNoSup(1,icond) = RefAmp;
        if doPlot,
            % waveforms
            figure(fh1);
            axes(hax(icond));
            xplot(xlim,RefAmp*[1 1],'g', 'linewidth',2);
            xplot(xlim,-RefAmp*[1 1],'g', 'linewidth',2);
            % Lissajous
            figure(fh2);
            axes(hax2(icond));
            xplot(xlim,RefAmp*[1 1],'g', 'linewidth',2);
        end
    end
    S.AmpShift = A2dB(abs(bsxfun(@rdivide,S.hrec,S.rmsNoSup)));
    S.PhaseShift = cunwrap(cangle(S.hrec./S.hrecNoSup));
    if doPlot && isfield(S, 'AmpShift'),
        fh3 = figure;
        set(fh3,'units', 'normalized', 'position', [0.148 0.259 0.589 0.578]);
        subplot(2,3,1);
        dplot(dt, S.AmpShift);
        ylabel('\Delta Amp (dB)');
        xlim([0,cycleDur]);
        title(Tstring, 'fontsize', 14, 'fontweight', 'bold', 'horizontalAlign', 'left');
        subplot(2,3,4);
        dplot(dt, S.PhaseShift);
        ylabel('\Delta Phase (cycle)');
        xlim([0,cycleDur]);
        xlabel('Time (ms)');
        subplot(2,3,[2 5]);
        plot(S.AmpShift, S.PhaseShift);
        xlabel('\Delta Amp (dB)');
        ylabel('\Delta Phase (cycle)');
        set(gca,'Xdir','reverse');
        subplot(2,3,3);
        plot(S.Fprobe/1e3, P2dB(mean(dB2P(S.AmpShift))), '-*', 'linewidth',2);
        xplot(S.Fprobe/1e3, min(S.AmpShift), '-');
        xplot(S.Fprobe/1e3, max(S.AmpShift), '-');
        ylabel('\Delta Amp (dB)');
        subplot(2,3,6);
        plot(S.Fprobe/1e3, mean(S.PhaseShift), '-*', 'linewidth',2);
        xplot(S.Fprobe/1e3, min(S.PhaseShift), '-');
        xplot(S.Fprobe/1e3, max(S.PhaseShift), '-');
        ylabel('\Delta Phase (dB)');
        xlabel('Fprobe (kHz)');
    end
end

%==========================================================
%==========================================================

function [lrec, hrec, cycleDur, dt] = local_Dsplit(D, irecChan, Fcutoff);
MaxNsidebands = 6;
[S, CFN, CP] = getcache([mfilename '_' expname(D)], {irec(D) MaxNsidebands});
if ~isempty(S),
    [lrec, hrec, cycleDur, dt] = deal(S.lrec, S.hrec, S.cycleDur, S.dt);
    return
end
GP = GenericStimparams(D);
St = D.Stim;
for icond = 1:GP.Ncond,
    [rec, dt] = anamean(D, irecChan, icond);
    [lrec(:,icond), hrec(:,icond), cycleDur] = local_split(rec, dt, icond, Fcutoff, GP, St, MaxNsidebands);
end
[Scal, Yunit] = conversionfactor(anachan(D, irecChan));
lrec = Scal*lrec;
hrec = Scal*hrec;
S = CollectInStruct(lrec, hrec, cycleDur, dt);
putcache(CFN, 500, CP, S, '-V6');

function [lrec, hrec, cycleDur] = local_split(rec, dt, icond, Fcutoff, GP, St, MaxNsidebands);
NsamOnset = round(GP.OnsetDelay(icond)/dt);
NsamRise = round(St.RiseDur/dt);
NsamFall = round(St.FallDur/dt);
NsamBurst = round(St.BurstDur/dt);
NsamCycle = St.NsamStimCycle;
% make sure that spectral analysis is commensurate with stimulus period
NsamTot = NsamCycle*floor(numel(rec)/NsamCycle); % largest integer # stim cycles fitting
rec = rec(1:NsamTot);
recDur = dt*NsamTot;
cycleDur = dt*NsamCycle;
df = 1e3/recDur; % freq spacing in Hz of spectrum
rec = simplegate(rec, NsamCycle);
Spec = fft(rec);
Nhalf = floor(numel(Spec)/2);
freq = Xaxis(Spec, df); % freq in Hz
Spec(freq>=St.Fsam/2) = 0; % positive freqs only
Spec = 1e6*i*Spec./freq; % BM velocity [mm/s]-> BM displacement [nm]
Spec(1) = 0;
% low freqs: suppressor & its harmonix; < cutoff
Lspec = Spec;
FSupHarmonix = (1:100)*St.Fsup_exact;
Lspec(~(local_freq_match(freq,FSupHarmonix))) = 0;
Lspec(freq>Fcutoff) = 0;
% high freqs: probe & sidebands > cutoff
Hspec = Spec;
Nsidebands = floor((St.Fprobe(icond)-Fcutoff)/St.Fsup_exact); % # sidebands to be included 
Nsidebands = min(Nsidebands, MaxNsidebands);
FprobeSidebands = St.Fprobe(icond) + (-Nsidebands:Nsidebands)*St.Fsup_exact;
Hspec(~(local_freq_match(freq,FprobeSidebands))) = 0;
Hspec(freq<Fcutoff) = 0;
lrec = 2*ifft(Lspec); % analytical waveform; LF
hrec = 2*ifft(Hspec); % analytical waveform; HF
NsamSteady = NsamBurst-NsamRise-NsamFall-2*NsamCycle;
NsamSteady = NsamCycle*floor(NsamSteady/NsamCycle);
isamOffset = NsamOnset+NsamRise+NsamCycle;
lrec = lrec(isamOffset+(1:NsamSteady));
hrec = hrec(isamOffset+(1:NsamSteady));
lrec = LoopMean(lrec, NsamCycle);
hrec = LoopMean(hrec, NsamCycle);
lrec = lrec-mean(lrec);

function qmatch = local_freq_match(allfreq, targetfreq);
% logical array indexing allfreq; true for ismember(allfreq, targetfreq)
% within rounding margins
allfreq = round(1e3*allfreq);
targetfreq = round(1e3*targetfreq);
qmatch = ismember(allfreq, targetfreq);

function local_labels(hax, Bh, Lh, D, St, DisplacementUnit, scaleFac);
for ah=Bh(:).',
    axes(ah);
    xlabel('Time (ms)');
end
for ah=Lh(:).',
    axes(ah);
    ylabel(['BM displacement (' DisplacementUnit ')']);
end

function local_Yscale(doScale, scaleFac, hax, Ncond, MAL);
if doScale,
    for icond = 1:Ncond,
        axes(hax(icond));
        ylim(1.1*scaleFac*MAL*[-1 1]);
    end
end

function [lrecRef, SPLRefLargeSup] = local_SupRef(D, RefLargeSup, irecChan, Fcutoff, Scal);
% use high-SPL suppr condition to get decent suppressor response
if isvoid(RefLargeSup),
    lrecRef = [];
    return;
end
[LrecRef, HrecRef, cycleDur, dt] = local_Dsplit(RefLargeSup, irecChan, Fcutoff);
StimD = D.Stim;
StimRef = RefLargeSup.Stim;
SPLRefLargeSup = StimRef.SupSPL;
dBdiff = StimD.SupSPL-SPLRefLargeSup;
icond = 1;
lrecRef = dB2A(dBdiff)*mean(LrecRef,2); % mean across conditions

function delay = local_delay(w1,w2,dt);
maxLag = round(numel(w1)/4);
w1 = repmat(w1(:),[5 1]);
w2 = repmat(w2(:),[5 1]);
[dum, ilag] = maxcorr(w2,w1,maxLag, 'abs');
delay = round(1e3*unique(ilag)*dt); % us



