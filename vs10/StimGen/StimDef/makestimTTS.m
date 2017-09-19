function P2=makestimTTS(P);
% MakestimTTS - stimulus generator for TTS stimGUI
%    P=MakestimTTS(P), where P is returned by GUIval, generates the stimulus
%    specified in P. MakestimTTS is typically called by StimGuiAction when
%    the user pushes the Check, Play or PlayRec button.
%    MakestimTTS does the following:
%        * Complete check of the stimulus parameters and their mutual
%          consistency, while reporting any errors
%        * Compute the stimulus waveforms
%        * Computation and broadcasting info about # conditions, total
%          stimulus duration, Max SPL, etc.
%
%    MakestimTTS renders P ready for D/A conversion by adding the following 
%    fields to P
%            Fsam: sample rate [Hz] of all waveforms. This value is
%                  determined by frequency content of the stimulus, but
%                  also by the Experiment definition P.Experiment, which
%                  may prescribe a minimum sample rate needed for ADC.
%          Fprobe: column array of probe freqs (Hz)
%      Fsup_exact: exact suppressor frequency (Hz) commensurate w sample rate
%   NsamStimCycle: number of samples in one stimulus cycle. Needed for data
%                  analysis using Fourier analysis.
%        Waveform: Waveform object array containing the samples in SeqPlay
%                  format.
%     Attenuation: scaling factors and analog attenuator settings for D/A
%    Presentation: struct containing detailed info on stimulus order,
%                  broadcasting of D/A progress, etc.
% 
%   See also stimdefFS.

P2 = []; % a premature return will result in []
if isempty(P), return; end
figh = P.handle.GUIfig;

% check & convert params. Note that helpers like evalfrequencyStepper
% report any problems to the GUI and return [] or false in case of problems.

% determine indiv probe freqs. Still approx values.
Fprobe = EvalFrequencyStepper(figh); 
if isempty(Fprobe), return; end
Ncond = size(Fprobe,1); % # conditions

% supp freq check
if P.SupFreq<P.Experiment.minStimFreq,
    GUImessage(figh, {'End frequency violates min stim frequency'...
        ['of ' num2str(P.Experiment.minStimFreq) ' Hz']},'error', 'SupFreq');
    return;
elseif P.SupFreq>P.Experiment.maxStimFreq,
    GUImessage(figh, {'Start frequency exceeds max stim frequency'...
        ['of ' num2str(P.Experiment.maxStimFreq) ' Hz']},'error', 'SupFreq');
    return;
end

% Durations & PlayTime; 
[okay, P]=EvalDurPanel(figh, P, Ncond);
if ~okay, return; end

% Sample rate & exact frequencies of probe & suppressor
P = local_exactfreqs(figh, P, Fprobe);
if isempty(P), 
    return;
end

% Generate the calibrated waveforms
P = local_wave(P); 

% Sort conditions, add baseline waveforms (!), provide info on varied parameter etc
P = sortConditions(P, 'Fprobe', 'Probe Frequency', 'Hz', P.StepFreqUnit);

% report durations
totBaseline = sum(SameSize(P.Baseline,[1 1])); % sum of pre- & post-stim baselines
Ttotal=ReportPlayTime(figh, Ncond, P.Nrep, P.ISI, totBaseline);

% Levels and active channels (must be called *after* adding the baseline waveforms)
[mxSPL P.Attenuation] = maxSPL(P.Waveform, P.Experiment);
okay=EvalSPLpanel(figh,P, mxSPL, P.Fprobe, '', {'SupSPL'});
if ~okay, return; end

P2=P;

%========================================

function P2 = local_exactfreqs(figh, P, Fprobe);
P2 = []; % allow premature return
doLock = isequal('to suppr.', P.LockProbe);
if doLock, % round Fprobe to nearest multiple of Fsup
    Frobe = P.SupFreq*round(Fprobe/P.SupFreq);
end
% sample rate needed for 5th order modulation in case of locking
P.Fsam = sampleRate(max([Fprobe(:); P.SupFreq])+doLock*5*P.SupFreq, P.Experiment); 
allokay = false; % pessimistic default
if isempty(P.Fsam) && doLock,
    GUImessage(figh, {'Highest probe freq and/or suppressor freq too high.', ...
        'Decrease highest probe freq, reduce suppr freq' ...
        'or give up locking the probe to the suppressor.'}, ...
        'error', {'StartFreq' 'EndFreq' 'SupFreq' 'LockProbe'});
elseif isempty(P.Fsam) && ~doLock,
    GUImessage(figh, {'Highest probe freq and/or suppressor freq too high.', ...
        'Decrease highest probe freq or reduce suppr freq.'}, 'error', ...
        {'StartFreq' 'EndFreq' 'SupFreq'});
elseif doLock && any(diff(sort(Frobe))==0),
    GUImessage(figh, {'Spacing of probe freqs is too narrow to', ...
        'realize locking to the suppressor.' ...
        'Decrease probe-frequency spacing, reduce suppr freq,' ...
        'or give up locking the probe to the suppressor.'}, ...
        'error', {'StepFreq' 'SupFreq' 'P.LockProbe'});
    return;
else, % no problems
    allokay = true; % passed all tests
end
if ~allokay, return; end
if doLock, % suppr cycle should contain integer # cycles
    P.Fsup_exact = P.Fsam/ceil(P.Fsam/P.SupFreq);
    P.Fprobe = P.Fsup_exact*round(Fprobe/P.Fsup_exact);
    P.NsamStimCycle = round(P.Fsam/P.Fsup_exact);
else, % round to nearest multiple of 1 Hz in TDT sense
    fiveTDTHz = P.Fsam/round(P.Fsam/5);
    P.Fsup_exact = fiveTDTHz*round(P.SupFreq/fiveTDTHz); 
    P.Fprobe = fiveTDTHz*round(Fprobe/fiveTDTHz); 
    P.NsamStimCycle = round(P.Fsam/fiveTDTHz);
end
P2 = P;


function P = local_wave(P);
Ncond = numel(P.Fprobe);
Chan = channelSelect(P.DAC, 'LR');
dt = 1e3/P.Fsam; % sample period in ms
P.StartPhase = [0; 0]; % start phases of probe and suppressor
for ichan = 1:numel(Chan),
    for icond=1:Ncond,
        freq = [P.Fprobe(icond); P.SupFreq]; % the 2 stim freqs
        [calibDL, calibDphi] = calibrate(P.Experiment, P.Fsam, Chan(ichan), freq);
        % waveform is generated @ the target SPL. Scaling is divided
        % between numerical scaling and analog attenuation later on.
        Amp = dB2A([P.SPL; P.SupSPL])*sqrt(2).*dB2A(calibDL); % numerical linear amplitudes of probe & suppr
        StoreDur = dt*P.NsamStimCycle; % dur of single stimulus cycle
        wtone = tonecomplex(Amp, freq, P.StartPhase, P.Fsam, StoreDur); % ungated waveform buffer; starting just after OnsetDelay
        if logical(P.Experiment.StoreComplexWaveforms); % if true, store complex analytic waveforms, real part otherwise
            wtone = wtone+ i*tonecomplex(Amp, freq, P.StartPhase+0.25, Fsam, StoreDur);
        end
        [Rise, Tail, NrepCyc] = cyclegate(wtone, dt, P.RiseDur, P.FallDur, P.BurstDur);
        NsamOnsetDelay = round(P.OnsetDelay/dt);
        P.Waveform(icond,ichan) = ...
            Waveform(P.Fsam, Chan(ichan), NaN, P.SPL, struct('fprobe', freq(1),'fsup', freq(2)), ...
            {0               Rise   wtone         Tail },...
            [NsamOnsetDelay  1      NrepCyc       1     ]);
        %   ^onset delay     ^rise  ^cyclic part  ^unfinished cycle+fall
        P.Waveform(icond,ichan) = AppendSilence(P.Waveform(icond,ichan), P.ISI); % pas zeros to ensure correct ISI
    end
end






