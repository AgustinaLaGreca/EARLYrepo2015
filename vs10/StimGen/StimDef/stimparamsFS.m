function P = stimparamsFS(S);
% stimparamsFS - generic stimulus parameters for FS stimulus protocol
%    G = stimparamsFS(S) takes the stimulus-defining struct S of an FS
%    protocol and returns the corresponding set of generic stimulus 
%    parameters im struct P. For a list of required fields in P, see
%    GenericStimparams. The trinity of functions stimdefFS, stimparamsFS,
%    and makestimFS completely define the FS stimulus protocol. 
%   
%    StimparamsFS is called by stimGUI during the preparation of the 
%    auditory stimuli (see Stimcheck). For older datasets, which do
%    not store the generic parameters, StimparamsFS may also be called by
%    GenericStimparams during data analysis.
%
%    See also GenericStimparams, stimGUI, stimDefPath, makestimFS, 
%    stimparamsFS.

Ncond = S.Presentation.Ncond;
dt = 1e3/S.Fsam; % sample period in ms

ID.StimType = S.StimType;
ID.Ncond = Ncond;
ID.Nrep  = S.Presentation.Nrep;
ID.Ntone = 1;
% ===backward compatibility
if ~isfield(S,'Baseline'), S.Baseline = 0; end
if ~isfield(S,'OnsetDelay'), S.OnsetDelay = 0; end
if ~isfield(S,'ModTheta'), S.ModTheta = 0; end
% ======timing======
T.PreBaselineDur = channelSelect('L', S.Baseline);
T.PostBaselineDur = channelSelect('R', S.Baseline);
T.ISI = ones(Ncond,1)*S.ISI;
T.BurstDur = channelSelect('B', S.Duration);
T.OnsetDelay = dt*floor(ones(Ncond,1)*S.OnsetDelay/dt); % always integer # samples
T.RiseDur = repmat(channelSelect('B', S.RiseDur), [Ncond 1]);
T.FallDur = repmat(channelSelect('B', S.FallDur), [Ncond 1]);
T.ITD = ones(Ncond,1)*S.ITD;
T.ITDtype = S.ITDtype;
T.TimeWarpFactor = ones(Ncond,1);
% ======freqs======
F.Fsam = S.Fsam;
F.Fcar = channelSelect('B', S.Fcar);
F.Fmod = channelSelect('B', S.Fmod);
F.LowCutoff = nan(Ncond,2);
F.HighCutoff = nan(Ncond,2);
F.FreqWarpFactor = ones(Ncond,1);
% ======startPhases & mod Depths
Y.CarStartPhase = repmat(channelSelect('B', S.WavePhase), [Ncond 1]);
Y.ModStartPhase = repmat(channelSelect('B', S.ModStartPhase), [Ncond 1]);
Y.ModDepth = repmat(channelSelect('B', S.ModDepth), [Ncond 1]);
Y.ModTheta = repmat(channelSelect('B', S.ModTheta), [Ncond 1]);
% ======levels======
L.SPL = repmat(channelSelect('B', S.SPL), [Ncond 1]);
L.SPLtype = 'per tone';
L.DAC = S.DAC;

P = structJoin(ID, '-Timing', T, '-Frequencies', F, '-Phases_Depth', Y, '-Levels', L);






