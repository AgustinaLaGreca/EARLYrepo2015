function P = toneStimDPOAE(P, varargin); 
% toneStim - compute tone stimulus
%   P = toneStim(P) computes the waveforms of tonal stimulus
%   The carrier frequencies are given by P.Fcar [Hz], which is a column 
%   vector (mono) or Nx2 vector (stereo). The remaining parameters are
%   taken from the parameter struct P returned by GUIval. The Experiment
%   field of P specifies the calibration, available sample rates etc
%   In addition to Experiment, the following fields of P are used to
%   compute the waveforms
%           Fcar: carrier frequency in Hz
%      WavePhase: starting phase (cycles) of carrier (0 means cos phase)
%    FreqTolMode: tolerance mode for freq rounding; equals exact|economic.
%        ModFreq: modulation frequency in Hz
%       ModDepth: modulation depth in %
%       ModStartPhase: modulation starting phase in cycles (0=cosine)
%       ModTheta: modulation angle in Cycle (0=AM, 0.25=QFM, other=mixed)
%            ISI: onset-to-onset inter-stimulus interval in ms
%     OnsetDelay: silent interval (ms) preceding onset (common to both DACs)
%       BurstDur: burst duration in ms including ramps
%        RiseDur: duration of onset ramp
%        FallDur: duration of offset ramp
%        FineITD: ITD imposed on fine structure in ms
%        GateITD: ITD imposed on gating in ms
%         ModITD: ITD imposed on modulation in ms
%            DAC: left|right|both active DAC channel(s)
%            SPL: carrier sound pressure level [dB SPL]
%
%   Most of these parameters may be a scalar or a [1 2] array, or 
%   a [Ncond x 1] or [Ncond x 2] or array, where Ncond is the number of 
%   stimulus conditions. The allowed number of columns (1 or 2) depends on
%   the character of the paremeter, i.e., whether it may have separate
%   values for the left and right DA channel. The exceptions are 
%   FineITD, GateITD, ModITD, and ISI, which which are allowed to have
%   only one column, and SPLtype and FreqTolMode, which are a single char 
%   strings that apply to all of the conditions.
%   
%   The output of toneStim is realized by updating/creating the following 
%   fields of P
%         Fsam:  sample rate [Hz] of all waveforms.
%         Fcar: adjusted to slightly rounded values to save memory using cyclic
%               storage (see CyclicStorage).
%         Fmod: modulation frequencies [Hz] in Ncond x Nchan matrix or column 
%               array. Might deviate slightly from user-specified values to
%               facilitate storage (see CyclicStorage).
%   CyclicStorage: the Ncond x Nchan struct with outputs of CyclicStorage
%     Duration: stimulus duration [ms] in Ncond x Nchan array. 
%    FineDelay: fine-structure-ITD realizing delays (columns denote Left,Right) 
%    GateDelay: gating-ITD realizing delays (columns denote Left,Right) 
%     ModDelay: modulation-ITD realizing delays (columns denote Left,Right) 
%     Waveform: Ncond x Nchan Waveform object containing the samples and 
%               additional info for D/A conversion.
%     GenericParamsCall: cell array for getting generic stimulus
%               parameters. Its value is 
%                   {@noiseStim struct([]) 'GenericStimParams'}
%               After substituting the updated stimulus struct for
%               struct([]), feval-ing this cell array will yield the 
%               generic stimulus parameters for noiseStim stimuli. 
%
%   For the realization of ITDs in terms of channelwise delays, see
%   ITD2delay.
%
%   toneStim(P, 'GenericStimParams') returns the generic stimulus
%   parameters for this class of tonal stimuli. This call is done via
%   GenericStimParams, based on the GenericParamsCall field described above.
%
%   See also makestimFS, SAMpanel, DurPanel, dePrefix, Waveform,
%   noiseStim, ITD2delay. 

if nargin>1 
    if isequal('GenericStimParams', varargin{1}),
        P = local_genericstimparams(P);
        return;
    else,
        error('Invalid second input argument.');
    end
end
S = [];
% test the channel restrictions described in the help text
error(local_test_singlechan(P,{'FineITD', 'GateITD', 'ModITD', 'ISI'}));
% There are Ncond=size(Fcar,1) conditions and Nch DA channels.
% Cast all numerical params in Ncond x Nch size, so we don't have to check
% sizes all the time.
[F1, F2, ModFreq, ModDepth, ModStartPhase, ModTheta, ...
    OnsetDelay, BurstDur, RiseDur, FallDur, ISI, WavePhase, ...
    FineITD, GateITD, ModITD, L1, L2, DAC] ...
    = SameSize(P.F1, P.F2, P.ModFreq, P.ModDepth, P.ModStartPhase, P.ModTheta, ...
    P.OnsetDelay, P.BurstDur, P.RiseDur, P.FallDur, P.ISI, P.WavePhase, ...
    P.FineITD, P.GateITD, P.ModITD, P.L1, P.L2, P.DAC);

% sign convention of ITD is specified by Experiment. Convert ITD to a nonnegative per-channel delay spec 
FineDelay = ITD2delay(FineITD(:,1), P.Experiment); % fine-structure binaural delay
GateDelay = ITD2delay(GateITD(:,1), P.Experiment); % gating binaural delay
ModDelay = ITD2delay(ModITD(:,1), P.Experiment); % modulation binaural delay

% find the single sample rate to realize all the waveforms while  ....
Fsam = sampleRate(F2, P.Experiment); % accounting for recording requirements minADCrate

% tolerances for memory-saving frequency roundings. 
[CarTol, ModTol] = FreqTolerance(F1, 0, P.FreqTolMode);

% compute # samples needed to store the waveforms w/o cyclic storage tricks
NsamTotLiteral = round(1e-3*sum(BurstDur)*Fsam);

NFreq = P.Ncond_XY(1);
NSPL = P.Ncond_XY(2);
% now compute the stimulus waveforms condition by condition, ear by ear.
% for ichan=1:Nchan,
%     chanStr = DAchanStr(ichan); % L|R

for ispl=1:NSPL
    
    for icond=1:NFreq,
        
        idx = icond + ((ispl-1)*NFreq);
        % evaluate cyclic storage to save samples
        C1 = CyclicStorage(F1(idx), 0, Fsam, BurstDur(idx), [CarTol(idx), ModTol(idx)], NsamTotLiteral(1));
        % compute the waveform
        
        [w1, fcar1] = local_Waveform('L', P.Experiment, Fsam, ISI(idx), ...
            FineDelay(idx), GateDelay(idx), ModDelay(idx), OnsetDelay(idx), RiseDur(idx), FallDur(idx), ...
            C1, WavePhase(idx), ModDepth(idx), ModStartPhase(idx), ModTheta(idx), L1(idx), P.StimType);
        w1 = setRep(w1,P.Nrep);
        
%         P.Waveform(idx,1) = w1;
%         % derived stim params
%         P.Fcar(idx,1) = fcar1;
%         P.CyclicStorage(idx,1) = C;
        
        
        C2 = CyclicStorage(F2(idx), 0, Fsam, BurstDur(idx), [CarTol(idx), ModTol(idx)], NsamTotLiteral(1));
        % compute the waveform
        
        [w2, fcar2] = local_Waveform('R', P.Experiment, Fsam, ISI(idx), ...
            FineDelay(idx), GateDelay(idx), ModDelay(idx), OnsetDelay(idx), RiseDur(idx), FallDur(idx), ...
            C2, WavePhase(idx), ModDepth(idx), ModStartPhase(idx), ModTheta(idx), L2(idx), P.StimType);
        w2 = setRep(w2,P.Nrep);
        
        
        switch DAC(1)
        case 'L'
            w1.Samples{1}= w1.Samples(:,1) + w2.Samples(:,1);
            w1 = AppendSilence(w1, ISI(idx));
            P.Waveform(idx,1) = w1;
            % derived stim params
            P.Fcar(idx,1) = fcar1;
            P.CyclicStorage(idx,1) = C1;
            
        case 'R'    %not implemented yet!!!!!!!!!!!        
            w2.Samples{1}= w1.Samples(:,1) + w2.Samples(:,1);
            w2 = AppendSilence(w2, ISI(idx));
            P.Waveform(idx,1) = w2;
            % derived stim params
            P.Fcar(idx,1) = fcar2;
            P.CyclicStorage(idx,1) = C2;
            
        otherwise          
            w1 = AppendSilence(w1, ISI(idx));
            w2 = AppendSilence(w2, ISI(idx));
            P.Waveform(idx,1) = w1;
            P.Waveform(idx,2) = w2;
            % derived stim params
            P.Fcar(idx,1) = fcar1;
            P.Fcar(idx,2) = fcar2;
            P.CyclicStorage(idx,1) = C1;
            P.CyclicStorage(idx,2) = C2;
        end
        
    end
end
P.Duration = SameSize(P.BurstDur, zeros(1,2)); 
P = structJoin(P, CollectInStruct(Fsam));
P.GenericParamsCall = {fhandle(mfilename) struct([]) 'GenericStimParams'};


%===================================================
%===================================================
function  [W, Fcar, Fmod] = local_Waveform(DAchan, EXP, Fsam, ISI, ...
    FineDelay, GateDelay, ModDelay, OnsetDelay, RiseDur, FallDur, ...
    C, WavePhase, ModDepth, ModStartPhase, ModTheta, SPL,StimType);
% Generate the waveform from the elementary parameters
%=======TIMING, DURATIONS & SAMPLE COUNTS=======
BurstDur = C.Dur;
% get sample counts of subsequent segments
SteadyDur = BurstDur-RiseDur-FallDur; % steady (non-windowd) part of tone
[NonsetDelay, NgateDelay, Nrise, Nsteady, Nfall] = ...
    NsamplesofChain([OnsetDelay, GateDelay, RiseDur, SteadyDur, FallDur], Fsam/1e3);
% For uniformity, cast literal storage in the form of a fake cyclic storage
useCyclicStorage = C.CyclesDoHelp && (Nsteady>=C.Nsam);
% if useCyclicStorage, % cyclic storage
%     NsamCyc = C.Nsam;  % # samples in cyclic buffer
%     NrepCyc = floor(Nsteady/NsamCyc); % # reps of cyclic buffer
%     NsamTail = rem(Nsteady,NsamCyc); % Tail containing remainder of cycles
%     Fcar = C.FcarProx; Fmod = C.FmodProx; % actual frequencies used in waveforms
% else, % literal storage: phrase as single rep of cyclic buffer + empty tail buffer
    NsamCyc = Nsteady;  % total # samples in steady-state portion
    NrepCyc = 1; % # reps of cyclic buffer
    NsamTail = 0; % No tail buffer
    Fcar = C.Fcar; Fmod = C.Fmod; % actual frequencies used in waveforms
% end
%=======FREQUENCIES, PHASES, AMPLITUDES and CALIBRATION=======
isModulated = ~isequal(0,ModDepth) && ~isequal(0,Fmod);
% SAM is implemented as 3-tone stim. Get the tone freqs.
if isModulated, freq = Fcar+[-1 0 1]*Fmod; % [Hz] lower sideband, carrier, upper sideband
else, freq = Fcar; % [Hz] only carrier
end
if strcmp(StimType,'SPONT') == 1
    calibDL = 0;
    calibDphi = 0;
else
    [calibDL, calibDphi] = calibrate(EXP, Fsam, DAchan, freq);
end
% waveform is generated @ the target SPL. Scaling is divided
% between numerical scaling and analog attenuation later on.
Amp = dB2A(SPL)*sqrt(2)*dB2A(calibDL); % numerical linear amplitudes of the carrier ...
if isModulated, Amp = Amp.*[ModDepth/200 1 ModDepth/200]; end %.. and sidebands if any
% Compute phase of numerical waveform at start of onset, i.e., first sample of rising portion of waveform.
StartPhase = WavePhase + calibDphi - 1e-3*FineDelay.*freq; % fine structure delay is realized in freq domain
if isModulated, % implement modulation phase by adjusting sideband phases
    % (1) Start phases [cycle] of the 3 tones at onset, excluding any ongoing delay
    StartPhase = StartPhase + [-1 0 1]*ModStartPhase; % 0/0.5 -> start at max/min envelope.
    % (2) Theta determines modulation type: 0=SAM; 0.25 = QFM; other values
    % produce mixed modulation (see vd Heijden & Joris, JARO 2010)
    StartPhase = StartPhase + [1 0 1]*ModTheta;
    % (3) whole waveform was already delayed a few lines above (using FineDelay). If the
    % requested modulation delay is different from the fine structure delay, we must
    % compensate for this.
    StartPhase = StartPhase - 1e-3*(ModDelay-FineDelay).*(freq-freq(2));
end

dt = 1e3/Fsam; % sample period in ms
% compute dur of stored tone, whether useCyclicStorage or not
StoreDur = dt*(NgateDelay + Nrise + (NsamCyc + NsamTail) + Nfall); % only a single cycled buf is used
wtone = tonecomplex(Amp, freq, StartPhase, Fsam, StoreDur); % ungated waveform buffer; starting just after OnsetDelay
if logical(EXP.StoreComplexWaveforms); % if true, store complex analytic waveforms, real part otherwise
    wtone = wtone+ i*tonecomplex(Amp, freq, StartPhase+0.25, Fsam, StoreDur);
end
wtone = ExactGate(wtone, Fsam, StoreDur-GateDelay, GateDelay, RiseDur, FallDur);
%set(gcf,'units', 'normalized', 'position', [0.519 0.189 0.438 0.41]); %xdplot(dt,wtone, 'color', rand(1,3)); error kjhkjh

% parameters stored w waveform. Mainly for debugging purposes.
NsamHead = NgateDelay + Nrise; % # samples in any gating delay + risetime portion
Nsam = CollectInStruct(NonsetDelay, NgateDelay, Nrise, Nsteady, Nfall, NsamHead);
Durs = CollectInStruct(BurstDur, RiseDur, FallDur, OnsetDelay, GateDelay, SteadyDur, FallDur);
Delays = CollectInStruct(FineDelay, GateDelay, ModDelay, OnsetDelay);
Param = CollectInStruct(C, Nsam, Durs, Delays, freq, Amp, StartPhase, SPL, useCyclicStorage);
% Patch together the segments of the tone, using the cycled storage format,
% or the fake version of it.
W = Waveform(Fsam, DAchan, NaN, SPL, Param, ...
    {0              wtone(1:NsamHead)  wtone(NsamHead+(1:NsamCyc))   wtone(NsamHead+NsamCyc+1:end)},...
    [NonsetDelay    1                  NrepCyc                       1], 1);
%    ^onset delay   ^gate_delay+rise   ^cyclic part                  ^remainder steady+fall  
 % pas zeros to ensure correct ISI
% AppendSlice was moved to the main function for this stimulus
%   The reason is that the samples needed to be added before the waveform was
%   appended. Similar method was used in DIZON.
 

function Mess = local_test_singlechan(P, FNS);
% test whether specified fields of P have single chan values
Mess = '';
for ii=1:numel(FNS),
    fn = FNS{ii};
    if size(P.(fn),2)~=1,
        Mess = ['The ''' fn ''' field of P struct must have a single column.'];
        return;
    end
end

function P = local_genericstimparams(S);
% extracting generic stimulus parameters. Note: this only works after
% SortCondition has been used to add a Presentation field to the
% stimulus-defining struct S.
Ncond = S.Presentation.Ncond;
dt = 1e3/S.Fsam; % sample period in ms
Nx1 = zeros(Ncond,1); % dummy for resizing
Nx2 = zeros(Ncond,2); % dummy for resizing
%
ID.StimType = S.StimType;
ID.Ncond = Ncond;
ID.Nrep  = S.Presentation.Nrep;
ID.Ntone = 1;
% ======timing======
T.PreBaselineDur = channelSelect('L', S.Baseline);
T.PostBaselineDur = channelSelect('R', S.Baseline);
T.ISI = SameSize(S.ISI, Nx1);
T.BurstDur = SameSize(channelSelect('B', S.Duration), Nx2);
T.OnsetDelay = SameSize(dt*floor(S.OnsetDelay/dt), Nx1); % always integer # samples
T.RiseDur = SameSize(channelSelect('B', S.RiseDur), Nx2);
T.FallDur = SameSize(channelSelect('B', S.FallDur), Nx2);
T.ITD = SameSize(S.ITD, Nx1);
T.ITDtype = S.ITDtype;
T.TimeWarpFactor = ones(Ncond,1);
% ======freqs======
F.Fsam = S.Fsam;
F.Fcar = SameSize(channelSelect('B', S.Fcar),Nx2);
F.Fmod = SameSize(channelSelect('B', S.ModFreq), Nx2);
F.F1 = SameSize(channelSelect('B', S.F1),Nx2);
F.F2 = SameSize(channelSelect('B', S.F2), Nx2);
F.LowCutoff = nan(Ncond,2);
F.HighCutoff = nan(Ncond,2);
F.FreqWarpFactor = ones(Ncond,1);
% ======startPhases & mod Depths
Y.CarStartPhase = nan([Ncond 2 ID.Ntone]);
Y.ModStartPhase = SameSize(channelSelect('B', S.ModStartPhase), Nx2);
Y.ModTheta = SameSize(channelSelect('B', S.ModTheta), Nx2);
Y.ModDepth = SameSize(channelSelect('B', S.ModDepth), Nx2);
% ======levels======
L.SPL = SameSize(channelSelect('B', S.SPL), Nx2);
L.L1 = SameSize(channelSelect('B', S.L1), Nx2);
L.L2 = SameSize(channelSelect('B', S.L2), Nx2);
L.SPLtype = 'per tone';
L.DAC = S.DAC;
P = structJoin(ID, '-Timing', T, '-Frequencies', F, '-Phases_Depth', Y, '-Levels', L);
P.CreatedBy = mfilename; % sign




