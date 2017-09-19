function P = stimparamsTTS(S);
% stimparamsTTS - generic stimulus parameters for TTS stimulus protocol
%    G = stimparamsTTS(S) takes the stimulus-defining struct S of an TTS
%    protocol and returns the corresponding set of generic stimulus 
%    parameters im struct P. For a list of required fields in P, see
%    GenericStimparams. The trinity of functions stimdefTTS, stimparamsTTS,
%    and makestimTTS completely define the TTS stimulus protocol. 
%   
%    StimparamsTTS is called by stimGUI during the preparation of the 
%    auditory stimuli (See Stimcheck). For older datasets, which do
%    not store the generic parameters, StimparamsTTS may also be called by
%    GenericStimparams during data analysis.
%
%    See also GenericStimparams, stimGUI, stimDefPath, makestimTTS, 
%    stimparamsTTS.

%public(S)
%P = 'waiver'; return;

Ncond = S.Presentation.Ncond;
dt = 1e3/S.Fsam; % sample period in ms

ID.StimType = S.StimType;
ID.Ncond = Ncond;
ID.Nrep  = S.Presentation.Nrep;
ID.Ntone = 2; % probe & suppressor

% ======timing======
T.PreBaselineDur = channelSelect('L', S.Baseline);
T.PostBaselineDur = channelSelect('R', S.Baseline);
T.ISI = ones(Ncond,1)*S.ISI;
T.BurstDur = repmat(channelSelect('B', S.BurstDur),[Ncond 1]);
T.OnsetDelay = dt*floor(ones(Ncond,1)*S.OnsetDelay/dt); % always integer # samples
T.RiseDur = repmat(channelSelect('B', S.RiseDur), [Ncond 1]); % ramp lasts 1 entire zwuis cycle
T.FallDur = T.RiseDur;
T.ITD = nan(Ncond,1);
T.ITDtype = 'none';
T.TimeWarpFactor = ones(Ncond,1);
% ======freqs======
F.Fsam = S.Fsam;
F.Fcar(:,1,:) = [S.Fprobe, repmat(S.Fsup_exact, [Ncond 1])]; F.Fcar = repmat(F.Fcar,[1 2 1]);
F.Fmod = zeros(Ncond,2);
F.LowCutoff = nan(Ncond,2);
F.HighCutoff = nan(Ncond,2);
F.FreqWarpFactor = ones(Ncond,1);
% ======startPhases & mod Depths
Y.CarStartPhase(1,1,:) = [S.StartPhase]; Y.CarStartPhase = repmat(Y.CarStartPhase, [Ncond 2 1]);
Y.ModStartPhase = zeros([Ncond 2]);
Y.ModTheta = zeros([Ncond 2]);
Y.ModDepth = zeros([Ncond 2]);
% ======levels======
L.SPL = repmat(channelSelect('B', S.SupSPL), [Ncond 1]); % suppressor level
L.SPLtype = 'per tone';
L.DAC = S.DAC;

P = structJoin(ID, '-Timing', T, '-Frequencies', F, '-Phases_Depth', Y, '-Levels', L);






