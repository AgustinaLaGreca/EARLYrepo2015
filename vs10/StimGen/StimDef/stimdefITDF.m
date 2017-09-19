function Params = stimdefITDF(EXP);
% stimdefITDF - definition of stimulus and GUI for ITDF stimulus paradigm
%    P=stimdefITDF(EXP) returns the definition for the ITDF (ITD&freq steps)
%    stimulus paradigm. The definition P is a GUIpiece that can be rendered
%    by GUIpiece/draw. Stimulus definition like stimmdefFS are usually
%    called by StimGUI, which combines the parameter panels with
%    a generic part of stimulus GUIs. The input argument EXP contains 
%    Experiment definition, which co-determines the realization of
%    the stimulus: availability of DAC channels, calibration, recording
%    side, etc.
%
%    See also stimGUI, stimDefDir, Experiment, makestimITDF.

PairStr = ' Pairs of numbers are interpreted as [left right].';
ClickStr = ' Click button to select ';
%==========Carrier frequency GUIpanel=====================
Fsweep = FrequencyStepper('carrier frequency', EXP);
ITDsweep = ITDstepper('ITD', EXP, '', 1); % last arg: specify different ITD types or not
% ---SAM
Sam = SAMpanel('SAM', EXP);

% ---Durations
Dur = DurPanel('-', EXP, '', 'basicsonly');
% ---Levels
Levels = SPLpanel('-', EXP);
% ---Pres
Pres = PresentationPanel_XY('Freq', 'ITD');
% ---Summary
summ = Summary(17);

%====================
Params=GUIpiece('Params'); % upper half of GUI: parameters

Params = add(Params, summ);
Params = add(Params, Fsweep, nextto(summ), [10 0]);
Params = add(Params, Sam, below(Fsweep), [0 10]);
Params = add(Params, ITDsweep, nextto(Fsweep), [10 0]);
Params = add(Params, Dur, below(ITDsweep) ,[-60 10]);
Params = add(Params, Levels, nextto(ITDsweep) ,[15 0]);
Params = add(Params, Pres, below(Levels) ,[35 25]);
Params = add(Params, PlayTime, below(Sam) , [30 5]);




