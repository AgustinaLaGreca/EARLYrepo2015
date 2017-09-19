function S = calibstim(Fsam, Flow, Fhigh, FMspeed, Amp);
% calibstim - generate FM sweep stimulus for calibration
%    S = calibstim(Fsam,Flow, Fhigh, FMspeed, Amp) returns calibration stimulus
%    consisting of a logarithmic FM sweep ranging from Flow Hz to Fhigh Hz,
%    sampled at Fsam Hz. FMspeed is the # octaves swept per second.
%    Amp is RMS of stimuli in Volt.

RampDur = 100; % ms ramps during which the stimulus is lingering @ end frequencies

Noct = log2(Fhigh/Flow); % total # octaves swept
SweepDur = round(1e3*Noct/FMspeed); % ms true sweep dur
NsamSweep = round(1e-3*SweepDur*Fsam); % # samples of sweep
NsamRamp =  round(1e-3*RampDur*Fsam); %  # samples of one ramp
InstFreq = logispace(Flow,Fhigh,NsamSweep).'; % instantaneous freq during sweep
InstFreq = [Flow*ones(NsamRamp-1,1); InstFreq; Fhigh*ones(NsamRamp,1)]; % include steady-state ramps
RadPerSam = 2*pi*InstFreq/Fsam; % convert Hz -> radians/sample
InstPhase = cumsum([0; RadPerSam]); % instantaneous phase in radians
Sweep = Amp*sqrt(2)*exp(i*InstPhase); % complex analytic waveform; RMS of real part is Amp
Sweep = localRamp(Sweep,round(NsamRamp/2)); % apply ramps

% --noise for estimate of overall delay---
SetRandState; addMLsig;
N = Gnoise(Flow, Fhigh, 0, 200, 0, [], 0, Fsam); % 200-ms white noise; same bandwidth as sweep
N = N.data; N = Amp*N/std(N); % normalize to RMS=Amp
N = localRamp(N,NsamRamp);

dt = 1e3/Fsam; % sample period in ms
Nsam = numel(Sweep);
% InstFeq now is freq between consecutive samples, not @ samples themselves. Extrapolate.
InstFreq = InstFreq([1 1:end end]);
InstFreq = 0.5*(InstFreq(1:end-1)+InstFreq(2:end));
S = CollectInStruct(Fsam, dt, Nsam, NsamRamp, Flow, Fhigh, FMspeed, Amp, SweepDur, RampDur, Sweep, InstFreq, N);


%======================
function W=localRamp(W,NsamRamp);
% supply cos^2 onset & offset ramps
Ramp = sin(linspace(0,pi/2,NsamRamp).');
W(1:NsamRamp) = Ramp.*W(1:NsamRamp);
W(end-NsamRamp+(1:NsamRamp)) = flipud(Ramp).*W(end-NsamRamp+(1:NsamRamp));



