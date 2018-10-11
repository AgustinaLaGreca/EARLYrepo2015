function P2=makestimNRHO(P);
% MakestimNRHO - stimulus generator for NRHO stimGUI
%    P=MakestimNRHO(P), where P is returned by GUIval, generates the stimulus
%    specified in P. MakestimNRHO is typically called by StimCheck when
%    the user pushes the Check, Play or PlayRec button.
%    MakestimNRHO does the following:
%        * Complete check of the stimulus parameters and their mutual
%          consistency, while reporting any errors
%        * Compute the stimulus waveforms
%        * Computation and broadcasting info about # conditions, total
%          stimulus duration, Max SPL, etc.
%
%    MakestimNRHO renders P ready for D/A conversion by adding the following 
%    fields to P
%            Fsam: sample rate [Hz] of all waveforms. This value is
%                  determined by the stimulus spectrum, but also by
%                  the Experiment definition P.Experiment, which may 
%                  prescribe a minimum sample rate needed for ADC.
%            Corr: column array of correaltion values realizing the steps.
%        Waveform: Waveform object array containing the samples in SeqPlay
%                  format.
%     Attenuation: scaling factors and analog attenuator settings for D/A
%    Presentation: struct containing detailed info on stimulus order,
%                  broadcasting of D/A progress, etc.
% 
%   See also stimdefNRHO.

P2 = []; % a premature return will result in []
if isempty(P), return; end
figh = P.handle.GUIfig;

% check & convert params. Note that helpers like evalphaseStepper
% report any problems to the GUI and return [] or false in case of problems.

% IPD: add it to stimparam struct P
P.Corr=EvalCorrStepper(figh, P); 
if isempty(P.Corr), return; end
Ncond = size(P.Corr,1); % # conditions

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
P.IPD = 0;
P = noiseStim(P); 

% Reduce storage of waveforms by setting some fields for play out
P.ReducedStorage = {'',''};
Chan = 'LR';
CorrChan = ipsicontra2LR(P.CorrChan, P.Experiment);
P.ReducedStorage{Chan~=CorrChan} = 'nonzero';
P.RX6_circuit = 'RX6seqplay-trig-2ADC';

% Sort conditions, add baseline waveforms (!), provide info on varied parameter etc
P = sortConditions(P, 'Corr', 'Interaural noise correlation', '', 'lin');

% Levels and active channels (must be called *after* adding the baseline waveforms)
[mxSPL, P.Attenuation] = maxSPL(P.Waveform, P.Experiment);
okay=EvalSPLpanel(figh,P, mxSPL, []);
if ~okay, return; end

% Summary
ReportSummary(figh, P);

P2=P;







