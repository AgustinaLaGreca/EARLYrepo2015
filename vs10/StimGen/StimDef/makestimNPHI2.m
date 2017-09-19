function P2=makestimNPHI(P);
% MakestimNPHI - stimulus generator for NPHI stimGUI
%    P=MakestimNPHI(P), where P is returned by GUIval, generates the stimulus
%    specified in P. MakestimNPHI is typically called by StimGuiAction when
%    the user pushes the Check, Play or PlayRec button.
%    MakestimFS does the following:
%        * Complete check of the stimulus parameters and their mutual
%          consistency, while reporting any errors
%        * Compute the stimulus waveforms
%        * Computation and broadcasting info about # conditions, total
%          stimulus duration, Max SPL, etc.
%
%    MakestimFS renders P ready for D/A conversion by adding the following 
%    fields to P
%            Fsam: sample rate [Hz] of all waveforms. This value is
%                  determined by the stimulus spectrum, but also by
%                  the Experiment definition P.Experiment, which may 
%                  prescribe a minimum sample rate needed for ADC.
%           Phase: column array of phases realizing the phase steps.
%        Waveform: Waveform object array containing the samples in SeqPlay
%                  format.
%     Attenuation: scaling factors and analog attenuator settings for D/A
%    Presentation: struct containing detailed info on stimulus order,
%                  broadcasting of D/A progress, etc.
%
%   Unlike the older version stimdefNPHI, stimdefNPHI2, uses the generic
%   waveform generator noisestim .
% 
%   See also stimdefNPHI2.

P2 = []; % a premature return will result in []
if isempty(P), return; end
figh = P.handle.GUIfig;

% check & convert params. Note that helpers like evalphaseStepper
% report any problems to the GUI and return [] or false in case of problems.

% IPD: add it to stimparam struct P
P.IPD=EvalPhaseStepper(figh, P); 
if isempty(P.IPD), return; end
Ncond = size(P.IPD,1); % # conditions

% split ITD in different types
[P.FineITD, P.GateITD, P.ModITD] = ITDparse(P.ITD, P.ITDtype);

% no heterodyning for this protocol
P.IFD = 0; % zero interaural frequency difference

% Noise parameters (SPL cannot be checked yet)
[okay, P.NoiseSeed] = EvalNoisePanel(figh, P);
if ~okay, return; end

% SAM (pass noise cutoffs to enable checking of out-of-freq-range sidebands)
okay=EvalSAMpanel(figh, P, [P.LowFreq P.HighFreq], {'LowFreqEdit' 'HighFreqEdit'});
if ~okay, return; end

% Durations & PlayTime messenger
okay=EvalDurPanel(figh, P, Ncond);
if ~okay, return; end

% Use generic noise generator to generate waveforms
P = noiseStim(P); 

% Sort conditions, add baseline waveforms (!), provide info on varied parameter etc
P = sortConditions(P, 'IPD', 'Interaural phase difference', 'cycle', 'lin');

% Levels and active channels (must be called *after* adding the baseline waveforms)
[mxSPL P.Attenuation] = maxSPL(P.Waveform, P.Experiment);
okay=EvalSPLpanel(figh,P, mxSPL, []);
if ~okay, return; end


P2=P;

%========================================
function P = local_wave(P);
NsamZeroBlock = 1e4; % # samples in silent block
ZeroBlock = zeros(NsamZeroBlock,1);
Ncond = size(P.Phase,1);
P.Fcar = nan(size(P.Phase));
P.Fmod = SameSize(P.ModFreq, P.Phase);
P.Dur = SameSize(max(P.BurstDur), zeros(Ncond,1)); % single column, see ToneStim
[P.Fsam, dt] = sampleRate(P.HighFreq+P.ModFreq, P.Experiment); % Fsam in Hz; dt in ms
binDelays = ITDsplit(P.ITD, P.ITDtype, P.Experiment, P.Experiment); % ITD converted into monaural delays
MinDur = P.BurstDur+ binDelays.maxGatingDelay; % burst dur enlarged w "worst case" additions from binaural calib delays & requested ITDs
binModDelays = ITDsplit(P.ModITD, 'ongoing', P.Experiment, P.Experiment); % convert AM ITDs to monaural delays. No effect on gating.
[AllChan, chanNum] = channelSelect(P.DAC, ['L' 'R'], 1:2);
[ModFreq, ModDepth, ModStartPhase, ModITD, ModTheta, LowFreq, HighFreq, SPL, BurstDur, RiseDur, FallDur] = ...
    SameSize(P.ModFreq, P.ModDepth, P.ModStartPhase, P.ModITD, P.ModTheta, P.LowFreq, P.HighFreq, ...
    P.SPL, P.BurstDur, P.RiseDur, P.FallDur, AllChan);
totNsam = round(P.ISI/dt);
for ichan=1:numel(AllChan),
    chanStr = AllChan(ichan); % L|R
    isRhoVarChannel = isequal(ear2DAchan(P.CorrUnit, P.Experiment),chanStr);
    % generate noise spectrum at "numerically correct" intensity
    if isRhoVarChannel, % realize interaural correlation by mixing
        NS = NoiseSpec(P.Fsam, MinDur, P.NoiseSeed, [LowFreq(ichan) HighFreq(ichan)], SPL(ichan), P.SPLUnit, P.Corr);
    else, % reference channel
        NS = NoiseSpec(P.Fsam, MinDur, P.NoiseSeed, [LowFreq(ichan) HighFreq(ichan)], SPL(ichan), P.SPLUnit, 1);
    end
    % apply ongoing delay & calibration while still in freq domain
    n = NS.Buf.*calibrate(P.Experiment, P.Fsam, chanStr, -NS.Freq, 1); % last one: complex phase factor; neg freqs: don't bother about freqs outside calib range
    odelay = binDelays.(chanStr).OngoingDelay; % ongoing delay (ms) for this DA chan L|R
    n = n.*exp(-2*pi*i*NS.Freq*1e-3*odelay); % apply this delay
    % go to time domain and apply AM & gating.
    n = ifft(n);
    n = n(1:ceil(MinDur/dt)); % throw away unused buffer tail
    if ModFreq(ichan)>0, 
        % apply AM. This is considered part of the "ongoing waveform". This means that the ongoing part of the binaural 
        % delays (see above) also applies to the AM. An exception is made when an explicit AM-ITD is specified;
        % in that case the ongoing delays from P.ITD  are overruled by the binaural AM delays.
        if isequal(0,P.ModITD), % use bin ITDS derived from P.ITD
            ph0 = ModStartPhase(ichan) - 1e-3*ModFreq(ichan)*binDelays.(chanStr).OngoingDelay;
        else, % use bin ITDS derived from P.ModITD
            ph0 = ModStartPhase(ichan) - 1e-3*ModFreq(ichan)*binModDelays.(chanStr).OngoingDelay;
        end
        n = SinMod(n, P.Fsam, ph0, ModFreq(ichan), ModDepth(ichan), ModTheta(ichan));
    end
    n = ExactGate(n,P.Fsam, BurstDur(ichan), binDelays.(chanStr).GatingDelay, RiseDur(ichan), FallDur(ichan));
    NheadZero = floor(P.OnsetDelay/dt);
    NtrailZero = totNsam-numel(n)-NheadZero;
    NremZero = rem(NtrailZero, NsamZeroBlock); 
    NrepZeroBlock = floor(NtrailZero/NsamZeroBlock);
    isVariedChannel = isequal(ear2DAchan(P.VariedChan, P.Experiment),chanStr);
    n = [n; zeros(NremZero,1)];
    for icond=1:Ncond,
        if isVariedChannel, 
            phasor = exp(2*pi*i*P.Phase(icond));
        else, phasor = 1;
        end
        P.Waveform(icond, ichan) = Waveform(P.Fsam, chanStr, max(abs(n)), SPL(ichan), ...
            CollectInStruct(binDelays, binModDelays), {0 real(phasor*n) ZeroBlock}, [NheadZero 1 NrepZeroBlock]);
    end
end

% # buffer blocks; 
% # comp & their freqs (also stop band)
% SPL of each cmp;
% set random seed
% compute Gaussian spectrum, 0 dBV, random phase (realize any ongoing ITD here!)
% get calibration data and make scalors to get correct
%   Ampl & phase of each cmp
% apply amp & phase factors
% ifft whole buffer
% apply SAM
% apply gating/ramps (realize any gating ITD here)
% apply phase factors
% take real part
% 








