function P2=makestimSup(P);
% MakestimSup - stimulus generator for SUP stimGUI
%    P=MakestimSUP(P), where P is returned by GUIval, generates the stimulus
%    specified in P. MakestimSup is typically called by StimGuiAction when
%    the user pushes the Check, Play or PlayRec button.
%    MakestimSup does the following:
%        * Complete check of the stimulus parameters and their mutual
%          consistency, while reporting any errors
%        * Compute the stimulus waveforms
%        * Computation and broadcasting info about # conditions, total
%          stimulus duration, Max SPL, etc.
%
%    MakestimSup renders P ready for D/A conversion by adding the following 
%    fields to P
%            Fsam: sample rate [Hz] of all waveforms. This value is
%                  determined by frequency content of the stimulus, but
%                  also by the Experiment definition P.Experiment, which
%                  may prescribe a minimum sample rate needed for ADC.
%          Fprobe: column array of probe freqs (Hz)
%             Sup: struct containing suppressor prescription. The
%                  suppressor is a Zwuis complex
%        Waveform: Waveform object array containing the samples in SeqPlay
%                  format.
%     Attenuation: scaling factors and analog attuater settings for D/A
%    Presentation: struct containing detailed info on stimulus order,
%                  broadcasting of D/A progress, etc.
% 
%   See also stimdefFS.

P2 = []; % a premature return will result in []
if isempty(P), return; end
figh = P.handle.GUIfig;

% check & convert params. Note that helpers like local_evalWarpStepper
% report any problems to the GUI and return [] or false in case of problems.

% extra duration check
if  P.SupDur > P.ISI,
    GUImessage(figh, {'Suppressor/probe duration exceeds ISI.'}, ...
        'error', {'ISI', 'SupDur'});
    return;
end

% determine indiv probe freqs. Still approx values.
Fprobe = EvalFrequencyStepper(figh); 
if isempty(Fprobe), return; end
Ncond = size(Fprobe,1); % # conditions

%---Supressor tone complex, preliminary version
Sup = local_evalsup(figh, P.Experiment, P);
if isempty(Sup),
    return;
end

% Exact frequencies
P = local_exactfreqs(figh, P, Sup, Fprobe);
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
okay=EvalSPLpanel(figh,P, mxSPL, P.Fprobe);
if ~okay, return; end


P2=P;

%========================================
function Sup = local_evalsup(figh, EXP, P, Fprobe);
%                   SupCenterFreq: 1500
%                      SupQuality: 'approx'
%                    SupBandWidth: 500
%                     SupBaseFreq: 2
%                          SupSPL: 65
%                          SupDur: 20000
%                        SupNcomp: 5
%                         SupSeed: 371786361
Sup = []; % allow premature returns
Nord = 3; % up to 3rd order
report(GUImessenger(figh,'SuppResult'), '');
if P.SupNcomp>1,
    DF = round(P.SupBandWidth/(P.SupNcomp-1));
else, DF = 0;
end
Tol = 0.2; % 20% freq tolerance
MaxNit = 100;
SetRandState(P.SupSeed);
try,
    [Fcmp, Nmult, DP, NtotIt] = ...
        uberzwuis(P.SupNcomp, Nord, P.SupBaseFreq, P.SupCenterFreq, DF, Tol, MaxNit);
catch,
    GUImessage(figh, lasterr, 'error');
    return;
end
report(GUImessenger(figh,'SuppResult'), [num2str(Nmult) ' overlaps']);
% {'SupCenterFreq' 'SupBandWidth'  'SupBaseFreq' 'SupNcomp' 'SupDur'};
allokay=0;
if isequal('perfect', P.SupQuality) && Nmult>0,
    GUImessage(figh, {'Perfect zwuis impossible with current param values.', ...
        'Relax the conditions or give up perfection.'}, 'error', ...
        {'SupBandWidth'  'SupBaseFreq' 'SupNcomp'});
    report(GUImessenger(figh,'SuppResult'), [num2str(Nmult) ' overlaps'], [0.8 0 0.1]);
elseif any(Fcmp<EXP.MinStimFreq),
    GUImessage(figh, {'Lowest suppressor component is too low.', ...
        'Increase center freq or reduce bandwidth.'}, 'error', ...
        {'SupCenterFreq' 'SupBandWidth'});
elseif any(Fcmp>EXP.MaxStimFreq),
    GUImessage(figh, {'Highest suppressor component is too high.', ...
        'Decrease center freq or reduce bandwidth.'}, 'error', ...
        {'SupCenterFreq' 'SupBandWidth'});
elseif P.SupDur<5e3/P.SupBaseFreq,
    GUImessage(figh, 'Duration must be at least 5 times base cycle.', 'error', ...
        {'SupDur' 'SupBaseFreq'});
else, allokay=1;
end
if ~allokay, return; end
Sup = CollectInStruct(Fcmp, Nmult, DP, NtotIt, '-', DF, Tol, MaxNit, Nord);

function P2 = local_exactfreqs(figh, P, Sup, Fprobe);
P2 = []; % allow premature return
ModFreqs = DiffMatrix(Sup.Fcmp); % beat freqs of suppressor
maxDPfreq = max(Fprobe) + max(ModFreqs(:)); % highest sideband created by suppression
minDPfreq = min(Fprobe) + min(ModFreqs(:)); % lowest sideband created by suppression
P.Fsam = sampleRate(maxDPfreq, P.Experiment); % sample rate needed for highest sideband
allokay = false; % pessimistic default
if isempty(P.Fsam),
    GUImessage(figh, {'Highest suppressor component is too high.', ...
        'Decrease center freq or reduce bandwidth.'}, 'error', ...
        {'SupCenterFreq' 'SupBandWidth'});
elseif minDPfreq<0,
    GUImessage(figh, {'Lowest suppressor component is negative.', ...
        'Increase center freq or reduce bandwidth.'}, 'error', ...
        {'SupCenterFreq' 'SupBandWidth'});
else, % no problems
    allokay = true; % passed all tests
end
if ~allokay, return; end
% now determine the exact freqencies for probe & zwuis stim. 
P.NsamZwuisPeriod = round(P.Fsam/P.SupBaseFreq); % integer # samples in single zwuis period
P.ZwuisBaseFreq = P.Fsam/P.NsamZwuisPeriod; % Hz exact zwuis base freq
P.Fprobe = P.ZwuisBaseFreq*round(Fprobe/P.ZwuisBaseFreq); % Hz exact probe freq: also integer # cycles
P.Fsup = P.ZwuisBaseFreq/P.SupBaseFreq*Sup.Fcmp; % correct zwuis freqs accordingly
P.Duration = P.SupDur;
P.SupDetails = Sup;
P2 = P;


function P = local_wave(P);
NsamZeroBlock = 1e4; % # samples in silent block
SetRandState(P.SupSeed);
NsamStimCycle = round(P.Fsam/P.ZwuisBaseFreq); 
riseWin = sin(linspace(0,pi/2,NsamStimCycle).').^2;
fallWin = flipud(riseWin);
Ncyc = floor(P.SupDur*1e-3*P.ZwuisBaseFreq); % # zwuis cycles in stimulus dur
% --suppressor waveform--
SPLperComp = P.SupSPL - P2dB(P.SupNcomp); % SPL of each sup component
P.SupStartPhase = rand(1,P.SupNcomp); % random start phases
[DL_sup, Dphi_sup, Dt, binCompDt] = calibrate(P.Experiment, P.Fsam, P.DAC, P.Fsup);
Time_s = (0:NsamStimCycle-1).'/P.Fsam; % time vector in s
SupWav = 0;
for icomp=1:P.SupNcomp,
    freq = P.Fsup(icomp);
    amp = sqrt(2)*dB2A(SPLperComp+DL_sup(icomp)); % linear amplitude
    ph0 = P.SupStartPhase(icomp)+Dphi_sup(icomp); % starting phase in cyles
    SupWav = SupWav + amp*cos(2*pi*(ph0+freq*Time_s));
end
% -- add probe waveforms--
[ichan, chanStr] = channelSelect(P.DAC,1:2, 'LR');
Ncond = numel(P.Fprobe);
DL_probe = calibrate(P.Experiment, P.Fsam, P.DAC, P.Fprobe);
NsamOnsetDelay = floor(1e-3*P.OnsetDelay*P.Fsam);
for icond = 1:Ncond,
    freq = P.Fprobe(icond);
    amp = sqrt(2)*dB2A(P.SPL+DL_probe(icond)); % linear amplitude
    probeWav = amp*cos(2*pi*freq*Time_s); % cos phase
    totWav = probeWav + SupWav;
    P.Waveform(icond, ichan) = Waveform(P.Fsam, chanStr, max(abs(totWav)), P.SPL(ichan), ...
        struct([]), {0 riseWin.*totWav totWav fallWin.*totWav }, [NsamOnsetDelay 1 Ncyc-2 1]); % 
    % append zeros to realize specified ISI
    P.Waveform(icond, ichan) = AppendSilence(P.Waveform(icond, ichan), P.ISI, NsamZeroBlock);
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






