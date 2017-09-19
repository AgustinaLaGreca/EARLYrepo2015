function Params = stimdefNPHI(EXP);
% stimdefNPHI - definition of stimulus and GUI for NPHI stimulus paradigm
%    P=stimdefNPHI(EXP) returns the definition for the NPHI (Noise with
%    variable phase) stimulus paradigm. The definition P is a GUIpiece that 
%    can be rendered by GUIpiece/draw. Stimulus definition like stimdefNPHI 
%    are usually called by StimGUI, which combines the parameter panels with
%    a generic part of stimulus GUIs. The input argument EXP contains 
%    Experiment definition, which co-determines the realization of
%    the stimulus: availability of DAC channels, calibration, recording
%    side, etc.
%
%    See also stimGUI, stimDefDir, Experiment, makestimNPHI.

PhiSweep = PhaseStepper('phase', EXP);
% ---SAM
Sam = SAMpanel('SAM', EXP, '', 1,1); % '': no prefix; 1,1: do include Theta and ITDs query
Sam = sameextent(Sam,PhiSweep,'X'); % adjust width of Mod to match PhiSweep
% ---Durations
Dur = DurPanel('-', EXP, '', 'nophase'); % exclude phase query in Dur panel
% ---Pres
Pres = PresentationPanel;
% Noise
Noise = NoisePanel('noise param', EXP);
%====================
Params=GUIpiece('Params'); % upper half of GUI: parameters
Params = add(Params, PhiSweep, 'below', [5 0]);
Params = add(Params, Sam, below(PhiSweep), [0 4]);
Params = add(Params, Noise, nextto(PhiSweep), [10 0]);
Params = add(Params, Dur, below(Noise) ,[0 20]);
Params = add(Params, Pres, nextto(Dur) ,[5 -20]);
Params = add(Params, PlayTime(), below(Dur) , [0 30]);




