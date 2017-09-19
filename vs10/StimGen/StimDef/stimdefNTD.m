function Params = stimdefNTD(EXP);
% stimdefNTD - definition of stimulus and GUI for NTD stimulus paradigm
%    P=stimdefNTD(EXP) returns the definition for the NTD (noise-ITD)
%    stimulus paradigm. The definition P is a GUIpiece that can be rendered
%    by GUIpiece/draw. Stimulus definition like stimmdefZW are usually
%    called by StimGUI, which combines the parameter panels with
%    a generic part of stimulus GUIs. The input argument EXP contains 
%    Experiment definition, which co-determines the realization of
%    the stimulus: availability of DAC channels, calibration, recording
%    side, etc.
%
%    See also stimGUI, stimDefDir, Experiment, makestimNTD.

% Noise
Noise = NoisePanel('noise param', EXP); % include SPL etc
ITD = ITDstepper('ITD', EXP, '', 1); % last arg: specify different ITD types or not
Sam = SAMpanel('modulation', EXP,'', 1, 0); % 1=include theta; 0 = exclude ITD
Dur = DurPanel('Durations', EXP, '', 'basicsonly');
Pres = PresentationPanel;
%====================
Params=GUIpiece('Params'); % upper half of GUI: parameters
Params = add(Params, Noise, 'below', [20 0]);
Params = add(Params, ITD, below(Noise));
Params = add(Params, Sam, nextto(ITD), [5 0]);
Params = add(Params, Dur, nextto(Noise), [20 0]);
Params = add(Params, Pres, below(Dur) ,[0 2]);
Params = add(Params, PlayTime(), below(ITD) ,[0 5]);




