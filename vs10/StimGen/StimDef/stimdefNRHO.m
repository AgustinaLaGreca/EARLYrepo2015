function Params = stimdefNRHO(EXP);
% stimdefNRHO - definition of stimulus and GUI for NRHO stimulus paradigm
%    P=stimdefNRHO(EXP) returns the definition for the NRHO (Noise with
%    variable phase) stimulus paradigm. The definition P is a GUIpiece that 
%    can be rendered by GUIpiece/draw. Stimulus definition like stimdefNRHO 
%    are usually called by StimGUI, which combines the parameter panels with
%    a generic part of stimulus GUIs. The input argument EXP contains 
%    Experiment definition, which co-determines the realization of
%    the stimulus: availability of DAC channels, calibration, recording
%    side, etc.
%
%    See also stimGUI, stimDefDir, Experiment, makestimNRHO.

CorrSweep = CorrStepper('noise correlation', EXP);
% ---SAM
Sam = SAMpanel('SAM', EXP, '', 1,0); % '': no prefix; 1,0: do include Theta, don't include ITDs query
CorrSweep = sameextent(CorrSweep,Sam, 'X'); % adjust width of Mod to match Sam
% ---Durations
Dur = DurPanel('-', EXP, '', 'nophase'); % exclude phase query in Dur panel
% ---Pres
Pres = PresentationPanel;
% Noise
Noise = NoisePanel('noise param', EXP, '', 'Corr');
% ---Summary
summ = Summary(17);

%====================
Params=GUIpiece('Params'); % upper half of GUI: parameters

Params = add(Params, summ);
Params = add(Params, CorrSweep, nextto(summ), [10 0]);
Params = add(Params, Sam, below(CorrSweep), [0 4]);
Params = add(Params, Noise, nextto(CorrSweep), [10 0]);
Params = add(Params, Dur, below(Noise) ,[0 20]);
Params = add(Params, Pres, nextto(Dur) ,[5 -20]);
Params = add(Params, PlayTime(), below(Dur) , [0 30]);




