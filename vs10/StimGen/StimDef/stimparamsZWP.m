function P = stimparamsZW(S);
% stimparamsFS - generic stimulus parameters for ZW stimulus protocol
%    G = stimparamsFS(S) takes the stimulus-defining struct S of an ZW
%    protocol and returns the corresponding set of generic stimulus 
%    parameters im struct P. For a list of required fields in P, see
%    GenericStimparams. The trinity of functions stimdefFS, stimparamsZW,
%    and makestimZW completely define the ZW stimulus protocol. 
%   
%    StimparamsZW is called by stimGUI during the preparation of the 
%    auditory stimuli (See Stimcheck). For older datasets, which do
%    not store the generic parameters, StimparamsZW may also be called by
%    GenericStimparams during data analysis.
%
%    See also GenericStimparams, stimGUI, stimDefPath, makestimZW, 
%    stimparamsZW.

%public(S)
Ncond = S.Presentation.Ncond;
dt = 1e3/S.Fsam; % sample period in ms

ID.StimType = S.StimType;
ID.Ncond = Ncond;
ID.Nrep  = S.Presentation.Nrep;
ID.Ntone = S.Ncomp;
% ===backward compatibility
if ~isfield(S,'Baseline'), S.Baseline = 0; end
if ~isfield(S,'OnsetDelay'), S.OnsetDelay = 0; end

% ======timing======
T.PreBaselineDur = channelSelect('L', S.Baseline);
T.PostBaselineDur = channelSelect('R', S.Baseline);
T.ISI = ones(Ncond,1)*S.ISI;
T.BurstDur = channelSelect('B', S.BurstDur);
T.OnsetDelay = dt*floor(ones(Ncond,1)*S.OnsetDelay/dt); % always integer # samples
T.RiseDur = repmat(channelSelect('B', S.RampDur), [Ncond 1]);
T.FallDur = repmat(channelSelect('B', S.RampDur), [Ncond 1]);
T.ITD = nan(Ncond,1);
T.ITDtype = 'none';
T.TimeWarpFactor = ones(Ncond,1);
% ======freqs======
F.Fsam = S.Fsam;
F.Fcar(1,1,:) = S.Fzwuis; F.Fcar = repmat(F.Fcar,[Ncond 2 1]);
F.Fmod = zeros(Ncond,2);
F.LowCutoff = nan(Ncond,2);
F.HighCutoff = nan(Ncond,2);
F.FreqWarpFactor = ones(Ncond,1);
% ======startPhases & mod Depths
Y.CarStartPhase(1,1,:) = S.StartPhase; Y.CarStartPhase = repmat(Y.CarStartPhase,[Ncond 2 1]);
Y.ModStartPhase = zeros([Ncond 2]);
Y.ModTheta = zeros([Ncond 2]);
Y.ModDepth = zeros([Ncond 2]);
% ======levels======
L.SPL = repmat(channelSelect('B', S.SPL), [Ncond 1]);
L.SPLtype = 'per tone';
L.DAC = S.DAC;

P = structJoin(ID, '-Timing', T, '-Frequencies', F, '-Phases_Depth', Y, '-Levels', L);






