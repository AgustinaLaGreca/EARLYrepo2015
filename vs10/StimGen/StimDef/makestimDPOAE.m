function P2=makestimDPOAE(P);
% MakestimFS - stimulus generator for FS stimGUI
%    P=MakestimFS(P), where P is returned by GUIval, generates the stimulus
%    specified in P. MakestimFS is typically called by StimGuiAction when
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
%                  determined by carrier & modulation freqs, but also by
%                  the Experiment definition P.Experiment, which may 
%                  prescribe a minimum sample rate needed for ADC.
%            Fcar: carrier frequencies [Hz] of all the presentations in an
%                  Nx2 matrix or column array
%        Waveform: Waveform object array containing the samples in SeqPlay
%                  format.
%     Attenuation: scaling factors and analog attuater settings for D/A
%    Presentation: struct containing detailed info on stimulus order,
%                  broadcasting of D/A progress, etc.
% 
%   See also toneStim, Waveform/maxSPL, Waveform/play, sortConditions, 
%   evalfrequencyStepper.

P2 = []; % a premature return will result in []
if isempty(P), return; end
figh = P.handle.GUIfig;

% Eval of Freq Stepper and finding F2 using the Ratio
x=EvalStepper(1,P.RatioSteps,P.RatioEnd);
P.F2=x*P.F1;
P.F2Unit='Hz';
Mess=['Frequency 2 (R): ' num2str(min(P.F2)) 'Hz to ' num2str(max(P.F2)) 'Hz'];
Mess=[Mess newline '                             in steps of ' num2str(min(P.F2)*P.RatioSteps) 'Hz'];
MM = GUImessenger(figh, ['F2Msg']);
report(MM,Mess);
if isempty(P.F1), return; end

L1=EvalStepper(P.L1Start,P.L1Step,P.L1End);
[P.F2, P.L1, P.Ncond_XY] = MixSweeps(P.F2, L1);
maxNcond = P.Experiment.maxNcond;
if prod(P.Ncond_XY)>maxNcond,
    Mess = {['Too many (>' num2str(maxNcond) ') stimulus conditions.'],...
        'Increase stepsize(s) or decrease range(s)'};
    GUImessage(figh, Mess, 'error', {'StartSPL' 'StepSPL' 'EndSPL' 'StartPhase' 'StepPhase' 'EndPhase'});
    return;
end
P.L1Unit='dB SPL';

P.L2=P.L1+P.DeltaL;
Mess=['L2          : from ' num2str(min(P.L2)) ' dB SPL'];
Mess=[Mess newline '                to '  num2str(max(P.L2)) ' dB SPL'];
MM = GUImessenger(figh, ['L2']);
report(MM,Mess);
if isempty(P.F1), return; end

% Durations & PlayTime; this also parses ITD/ITDtype and adds ...
P.ITD=0;
[okay, P]=EvalDurPanel(figh, P, P.Ncond_XY);% ... FineITD, GateITD, ModITD fields to P
if ~okay, return; end

% Different logic is used and variable names are changed, dealing 0 to
% bypass stim parameters check
[P.Fcar, P.ModFreq, P.SPL] = deal(0);

% No modulation panel
[P.ModFreq, P.ModDepth, P.ModStartPhase, P.ModITD, P.ModTheta, P.GateDelay, P.WavePhase]=deal(0);

% This stimuslus does not have ITD
[P.FineITD, P.GateITD, P.ModITD] = deal(0);

% Determine sample rate and actually generate the calibrated waveforms
P = toneStimDPOAE(P); % P contains both Experiment (calib etc) & params, including P.Fcar 

% Process visiting order of stimulus conditions
VisitOrder = EvalPresentationPanel_XY(figh, P, P.Ncond_XY);
if isempty(VisitOrder), return; end

% Sort conditions, add baseline waveforms (!), provide info on varied parameter etc
P.DAC='Both';
P = sortConditions(P, {'F2','L1'}, {'Freq 2','SPL 1'}, {'Hz','dB'}, {P.F2Unit,P.L1Unit});
% P = sortConditions(P, 'F1', 'Freq 1','Hz', P.F1Unit);

% Levels and active channels (must be called *after* adding the baseline waveforms)
[mxSPL P.Attenuation] = maxSPL(P.Waveform, P.Experiment);
okay=EvalSPLpanelDPOAE(figh,P, mxSPL, P.Fcar);
if ~okay, return; end
% Summary
ReportSummaryDPOAE(figh, P);

% 'TESTING MAKEDTIMFS'
% P.Duration
% P.Duration = []; % 
% P.Fcar = [];

% everything okay: return P
P2=P;











