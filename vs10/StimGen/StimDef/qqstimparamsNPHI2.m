function P = stimparamsNPHI(S);
% stimparamsNPHI - generic stimulus parameters for NPHI stimulus protocol
%    G = stimparamsNPHI(S) takes the stimulus-defining struct S of an NPHI
%    protocol and returns the corresponding set of generic stimulus 
%    parameters im struct P. For a list of required fields in P, see
%    GenericStimparams. The trinity of functions stimdefNPHI, stimparamsNPHI,
%    and makestimNPHI completely define the NPHI stimulus protocol. 
%   
%    StimparamsNPHI is called by stimGUI during the preparation of the 
%    auditory stimuli (see Stimcheck). For older datasets, which do
%    not store the generic parameters, StimparamsNPHI may also be called by
%    GenericStimparams during data analysis.
%
%    See also GenericStimparams, stimGUI, stimDefPath, makestimNPHI, 
%    stimparamsNPHI.

%public(S)
Ncond = S.Presentation.Ncond;
dt = 1e3/S.Fsam; % sample period in ms

ID.StimType = S.StimType;
ID.Ncond = Ncond;
ID.Nrep  = S.Presentation.Nrep;
ID.Ntone = 0;
% ===backward compatibility
if ~isfield(S,'Baseline'), S.Baseline = 0; end
if ~isfield(S,'OnsetDelay'), S.OnsetDelay = 0; end

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
F.Fcar = nan(Ncond,2,ID.Ntone);
F.Fmod = SameSize(channelSelect('B', S.ModFreq), zeros(Ncond,1));
F.LowCutoff = repmat(channelSelect('B', S.LowFreq), [Ncond 1]);
F.HighCutoff = repmat(channelSelect('B', S.HighFreq), [Ncond 1]);
F.FreqWarpFactor = ones(Ncond,1);
% ======startPhases & mod Depths
Y.CarStartPhase = nan([Ncond 2 ID.Ntone]);
Y.ModStartPhase = repmat(channelSelect('B', S.ModStartPhase), [Ncond 1]);
Y.ModDepth = repmat(channelSelect('B', S.ModDepth), [Ncond 1]);
% ======levels======
L.SPL = repmat(channelSelect('B', S.SPL), [Ncond 1]);
if isequal('dB/Hz', S.SPLUnit), L.SPLtype = 'spectrum level';
elseif isequal('dB SPL', S.SPLUnit), L.SPLtype = 'total level';
end
L.DAC = S.DAC;

P = structJoin(ID, '-Timing', T, '-Frequencies', F, '-Phases_Depth', Y, '-Levels', L);






