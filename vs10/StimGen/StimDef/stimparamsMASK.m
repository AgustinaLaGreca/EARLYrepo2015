function P = stimparamsMASK(S);
% stimparamsFS - generic stimulus parameters for MASK stimulus protocol
%    G = stimparamsFS(S) takes the stimulus-defining struct S of an MASK
%    protocol and returns the corresponding set of generic stimulus 
%    parameters im struct P. For a list of required fields in P, see
%    GenericStimparams. The trinity of functions stimdefFS, stimparamsMASK,
%    and makestimMASK completely define the MASK stimulus protocol. 
%   
%    StimparamsMASK is called by stimGUI during the preparation of the 
%    auditory stimuli (See Stimcheck). For older datasets, which do
%    not store the generic parameters, StimparamsMASK may also be called by
%    GenericStimparams during data analysis.
%
%    Note 1: MASK has two simultaneous stimulus components (Noise & tone). 
%    The choice of most parameters has been made in favor of the noise
%    waveform. See code for details.
%
%    Note 2: MASK applies time warping. The generic paramers take warping
%    into account, e.g., the ISI varies across conditions. This is why the
%    generic stimulus parameters look different from the values enetered in
%    the stimulus GUI, which are "pre-warp."
%
%    See also GenericStimparams, stimGUI, stimDefPath, makestimMASK, 
%    stimparamsMASK.

%public(S)
Ncond = S.Presentation.Ncond;
dt = 1e3/S.Fsam; % sample period in ms

ID.StimType = S.StimType;
ID.Ncond = Ncond;
ID.Nrep  = S.Presentation.Nrep;
ID.Ntone = 1;
% ===backward compatibility
if ~isfield(S,'Baseline'), S.Baseline = 0; end
if ~isfield(S,'MnoiseOnsetDelay'), S.MnoiseOnsetDelay = 0; end

% ======timing======
Twarp = 2.^(-S.Warp); % time scaling
T.PreBaselineDur = channelSelect('L', S.Baseline);
T.PostBaselineDur = channelSelect('R', S.Baseline);
T.ISI = Twarp*S.ISI;
T.BurstDur = channelSelect('B', S.Duration);
T.OnsetDelay = dt*floor(Twarp*S.MnoiseOnsetDelay/dt); % always integer # samples
T.RiseDur = Twarp*channelSelect('B', S.MnoiseRiseDur);
T.FallDur = Twarp*channelSelect('B', S.MnoiseFallDur);
T.ITD = repmat(S.MnoiseITD, [Ncond 1]);
T.ITDtype = S.MnoiseITDtype;
T.TimeWarpFactor = Twarp;
% ======freqs======
Fwarp = 1./Twarp; % spectral scaling factor
F.Fsam = S.Fsam;
F.Fcar = channelSelect('B', S.Fcar);
F.Fmod = zeros(Ncond,2);
F.LowCutoff = Fwarp*channelSelect('B', S.MnoiseLowFreq);
F.HighCutoff = Fwarp*channelSelect('B', S.MnoiseHighFreq);
F.FreqWarpFactor = Fwarp;
% ======startPhases & mod Depths
Y.CarStartPhase = zeros([Ncond 2]);
Y.ModStartPhase = zeros([Ncond 2]);
Y.ModTheta = zeros([Ncond 2]);
Y.ModDepth = zeros([Ncond 2]);
% ======levels======
L.SPL = repmat(channelSelect('B', S.MnoiseSPL), [Ncond 1]);
if isequal('dB/Hz', S.MnoiseSPLUnit), L.SPLtype = 'spectrum level';
elseif isequal('dB SPL', S.MnoiseSPLUnit), L.SPLtype = 'total level';
end
L.DAC = S.DAC;

P = structJoin(ID, '-Timing', T, '-Frequencies', F, '-Phases_Depth', Y, '-Levels', L);






