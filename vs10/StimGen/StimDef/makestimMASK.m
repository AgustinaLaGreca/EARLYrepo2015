function P2=makestimMask(P);
% MakestimMASK - stimulus generator for MASK stimGUI
%    P=MakestimMASK(P), where P is returned by GUIval, generates the stimulus
%    specified in P. MakestimNPHI is typically called by StimGuiAction when
%    the user pushes the Check, Play or PlayRec button.
%    MakestimFS does the following:
%        * Complete check of the stimulus parameters and their mutual
%          consistency, while reporting any errors
%        * Compute the stimulus waveforms
%        * Computation and broadcasting info about # conditions, total
%          stimulus duration, Max SPL, etc.
%
%    MakestimMASK renders P ready for D/A conversion by adding the following 
%    fields to P
%            Fsam: sample rate [Hz] of all waveforms. This value is
%                  determined by frequency content of the stimulus, but
%                  also by the Experiment definition P.Experiment, which
%                  may prescribe a minimum sample rate needed for ADC.
%            Warp: column array of warp values realizing the warp steps.
%        Waveform: Waveform object array containing the samples in SeqPlay
%                  format.
%     Attenuation: scaling factors and analog attuater settings for D/A
%    Presentation: struct containing detailed info on stimulus order,
%                  broadcasting of D/A progress, etc.
% 
%   See also stimdefNPHI.

P2 = []; % a premature return will result in []
if isempty(P), return; end
figh = P.handle.GUIfig;

% check & convert params. Note that helpers like local_evalWarpStepper
% report any problems to the GUI and return [] or false in case of problems.

% determine indiv warp values & add them to stimparam struct P
[P.Warp, Ncond, Mess] = local_evalWarpStepper(figh, P, P.Experiment); 
if isempty(P.Warp), return; end

% Noise parameters (SPL cannot be checked yet)
okay = EvalNoisePanel(figh, P, 'Mnoise');
if ~okay, return; end

% Durations & PlayTime messenger
WarpFactor = 2.^P.Warp; % spectral warp factor: octaves -> freq ratio
meanISI = P.ISI*mean(1./WarpFactor); % duration warp factor is reciprocal of spectral warp factor
okay=EvalDurPanel(figh, P, Ncond, 'Mnoise', meanISI); 
if ~okay, return; end

% extra duration check
BurstStart = min(P.ToneOnset,0); % noise always starts at t=0
BurstEnd = max(P.ToneOnset + P.ToneDur, P.MnoiseBurstDur+abs(P.MnoiseITD));
BurstDur = BurstEnd-BurstStart;
if any(BurstDur > P.ISI),
    GUImessage(figh, {'Duration exceeds ISI.', 'Note that noise ITD and tone onset affect the duration.'}, ...
        'error', {'ISI', 'ToneOnset', 'ToneDur', 'MnoiseBurstDur', 'MnoiseITD'});
    return;
end

% Determine sample rate and generate the calibrated waveforms
P = local_wave(figh,P); 
if isempty(P), return; end

% Sort conditions, add baseline waveforms (!), provide info on varied parameter etc
P = sortConditions(P, 'Warp', 'Warp value', 'octave', 'lin');

% Levels and active channels (must be called *after* adding the baseline waveforms)
[mxSPL P.Attenuation] = maxSPL(P.Waveform, P.Experiment);
okay=EvalSPLpanel(figh,P, mxSPL, P.Fcar, 'Mnoise');
if ~okay, return; end


P2=P;

%========================================
function [Warp, Ncond, Mess] = local_evalWarpStepper(figh, P, EXP); 
[Warp, Mess]=EvalStepper(P.StartWarp, P.StepWarp, P.EndWarp, 'linear', ...
    P.AdjustWarp, [-inf inf], EXP.maxNcond); % +/- inf: no a priori limits imposed to warp values
if isequal('nofit', Mess),
    Mess = {'Stepsize does not exactly fit Warp bounds', ...
        'Adjust Warp parameters or toggle Adjust button.'};
elseif isequal('largestep', Mess)
    Mess = 'Step size exceeds Warp range';
elseif isequal('cripple', Mess)
    Mess = 'Different # Warp steps in the two DA channels';
elseif isequal('toomany', Mess)
    Mess = {['Too many (>' num2str(EXP.maxNcond) ') Warp steps.'] 'Increase stepsize or decrease range'};
end
GUImessage(figh, Mess, 'error', {'StartWarp', 'StepWarp', 'EndWarp'});
Ncond = size(Warp,1);


function P = local_wave(figh,P);
NsamZeroBlock = 1e4; % # samples in silent block
UpsampleFactor = 2;
NsamWarpMargin = 1000; % # extra samples needed to avoid boundary effects when oversampling (see below)
WarpFactor = 2.^P.Warp;
[MinWarpFac, MaxWarpFac] = minmax(WarpFactor);
Ncond = size(P.Warp,1);
[AllChan, chanNum] = channelSelect(P.DAC, ['L' 'R'], 1:2);
Nchan = numel(AllChan);
Ftone = SameSize(P.ToneFreq, zeros(Ncond, Nchan)); % true tone freqs, accounting for warp
P.Fcar = WarpFactor*ones(1,Nchan).*Ftone;
P.Fmod = deal(zeros(size(P.Warp)));
P.Duration = max([P.MnoiseBurstDur (P.ToneOnset+P.ToneDur)])./WarpFactor; % single column, see ToneStim
P.DurOkay = 1;
[P.Fsam, dt] = sampleRate(MaxWarpFac*[P.MnoiseHighFreq P.ToneFreq], P.Experiment); % Fsam in Hz; dt in ms
if isempty(P.Fsam),
    Mess = {'Frequencies too high to be realized.', 'Note that warping affects the spectrum.'};
    GUImessage(figh, Mess, 'error', {'StartWarp', 'StepWarp', 'EndWarp' 'MnoiseHighFreq' 'ToneFreq'});
end
% derived params: unwarped noise first.
binDelays = ITDsplit(P.MnoiseITD, P.MnoiseITDtype, P.Experiment, P.Experiment); % ITD converted into monaural delays
MinBufDur = P.MnoiseBurstDur + binDelays.maxGatingDelay + NsamWarpMargin*dt; % burst dur enlarged w "worst case" additions from binaural calib delays & requested ITDs
[LowFreq, HighFreq, SPL, BurstDur, RiseDur, FallDur, NoisePhase] = ...
    SameSize(P.MnoiseLowFreq, P.MnoiseHighFreq, ...
    P.MnoiseSPL, P.MnoiseBurstDur, P.MnoiseRiseDur, P.MnoiseFallDur, P.MnoiseWavePhase, ...
    AllChan);
% derived params: unwarped tone next.
toneBinDelays = ITDsplit(0, 'waveform', P.Experiment, P.Experiment); % zero ITD converted into monaural delays
[ToneFreq, ToneSPL, ToneOnset, ToneDur, ToneRiseDur, ToneFallDur] ...
    = SameSize(P.ToneFreq, P.ToneSPL, P.ToneOnset, P.ToneDur, P.ToneRiseDur, P.ToneFallDur, AllChan);
% who starts first, noise or tone? The first one start at t=0 by def. Adapt
% params accordingingly.
if any(ToneOnset<0), % noise must be delayed re tone onset. 
    NoiseOnset = -min(ToneOnset)*ones(size(ToneOnset));
    MinBufDur = max(MinBufDur + NoiseOnset);
    ToneOnset = ToneOnset+NoiseOnset;
else, NoiseOnset = [0 0];
end
for ichan=1:numel(AllChan),
    chanStr = AllChan(ichan); % L|R
    isRhoChannel = isequal(ear2DAchan(P.MnoiseCorrUnit, P.Experiment),chanStr);
    % generate noise spectrum at "numerically correct" intensity
    if isRhoChannel, rho = P.MnoiseCorr; % realize interaural correlation by mixing
    else, rho=1; % reference channel
    end
    NS = NoiseSpec(P.Fsam, MinBufDur, P.MnoiseConstNoiseSeed, [LowFreq(ichan) HighFreq(ichan)], SPL(ichan), P.MnoiseSPLUnit, rho);
    % apply ongoing delay & calibration while still in freq domain
    n = NS.Buf.*calibrate(P.Experiment, P.Fsam, chanStr, -NS.Freq, 1); % last one: complex phase factor; neg freqs: totrate out-of range freqs
    odelay = NoiseOnset(ichan) + binDelays.(chanStr).OngoingDelay; % ongoing delay (ms) for this DA chan L|R
    n = n.*exp(-2*pi*i*NS.Freq*1e-3*odelay); % apply this delay
    % go to time domain, reduce length, apply phase & take real part
    n = ifft(n);
    n = n([end-NsamWarpMargin+1:end 1:ceil(MinBufDur/dt)+NsamWarpMargin]); % throw away unused buffer tail, but inclusing flanking margins
    Ph0 = NoisePhase(ichan);
    n = real(exp(2*pi*i*Ph0)*n);
    n = resample(n, UpsampleFactor, 1); % oversample in anticipation of warps
    dtn = 1e3/P.Fsam/UpsampleFactor; % new, smaller sample period
    n([1:UpsampleFactor*NsamWarpMargin (end+1-UpsampleFactor*NsamWarpMargin):end]) = []; % remove margins
    n = ExactGate(n,UpsampleFactor*P.Fsam, BurstDur(ichan), ...
        NoiseOnset(ichan) + binDelays.(chanStr).GatingDelay, RiseDur(ichan), FallDur(ichan));
    % tone
    ToneCalibFactor = calibrate(P.Experiment, P.Fsam, chanStr, ToneFreq(ichan), 1); % last one: complex phase factor
    TotToneOnset = ToneOnset(ichan)+toneBinDelays.(chanStr).Acoust_Comp; % onset delay w acoustic compensation included
    Tone = localTone(dtn, ToneCalibFactor, ToneFreq(ichan), NS.specSPL+ToneSPL(ichan), ...
        TotToneOnset, ToneDur(ichan), ToneRiseDur(ichan), ToneFallDur(ichan));
    n = padandadd(n,Tone); % trailing zeros are padded to compensate for any length differences
    % the conditions differ in their warping, realized by linear
    % interpolation, followed by downsampling.
    for icond=1:Ncond,
        wf = WarpFactor(icond); % spectral warp factor; wf>1 -> freq scaled up -> dur reduced
        isi = P.ISI/wf; % true ISI after warping
        dur = dtn*numel(n); % current dur of noise+tone
        nsam = numel(n);
        newdur = dur/wf; % warped duration of noise+tone
        newnsam = round(newdur/dtn);
        nw = interp1(linspace(0,1,nsam), n, linspace(0,1,newnsam).');
        nw = resample(nw,1,UpsampleFactor); % downsample toward original sample rate
        NsamOnsetDelay = floor(P.MnoiseOnsetDelay/dt/wf); % # zero samples to be prepended
        wparam = CollectInStruct(binDelays, toneBinDelays, isi); % params passed to waveform
        P.Waveform(icond, ichan) = Waveform(P.Fsam, chanStr, max(abs(nw)), SPL(ichan), ...
            wparam, {0 nw}, [NsamOnsetDelay 1]); % last 1 is Nrep
        %ddd= duration(P.Waveform(icond, ichan))*wf
        %isiwf = isi*wf
        %disp('========')
        P.Waveform(icond, ichan) = AppendSilence(P.Waveform(icond, ichan), isi, NsamZeroBlock); % add zeros
    end
end

function Tone = localTone(dt, CalibFactor, Freq, SPL, Onset, Dur, RiseDur, FallDur);
% Compute tone to be added to calibration-compensated noise waveform.
% cos phase at start of rise window. First sample of the tone is to be
% added to first sample of noise. This alignment determines the onset.
Nsam = ceil((Onset+Dur)/dt); % # samples of ungated tone waveform
ph0_rad = -2*pi*(Onset*1e-3*Freq); % starting phase in radians at start of n
omega = 2*pi*1e-3*Freq*dt; % tone frequency in radians per sample
Tone = dB2A(SPL)*exp(i*(ph0_rad+omega*(0:Nsam-1))).'; % complex tone, long enough to be gated
Tone = real(CalibFactor*Tone); % apply calibration fine structure & restrict to lin part
Tone = ExactGate(Tone, 1e3/dt, Dur, Onset, RiseDur, FallDur);






