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
%   See also stimdefNPHI.

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





