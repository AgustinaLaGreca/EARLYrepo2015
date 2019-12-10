function Params = stimdefDPOAE(EXP);
% stimdefFS - definition of stimulus and GUI for FS stimulus paradigm
%    P=stimdefFS(EXP) returns the definition for the FS (Freqency Sweep)
%    stimulus paradigm. The definition P is a GUIpiece that can be rendered
%    by GUIpiece/draw. Stimulus definition like stimmdefFS are usually
%    called by StimGUI, which combines the parameter panels with
%    a generic part of stimulus GUIs. The input argument EXP contains 
%    Experiment definition, which co-determines the realization of
%    the stimulus: availability of DAC channels, calibration, recording
%    side, etc.
%   
%    See also stimGUI, stimDefDir, Experiment, makestimFS, stimparamsFS.

PairStr = ' Pairs of numbers are interpreted as [left right].';
ClickStr = ' Click button to select ';
%==========Carrier frequency GUIpanel=====================
FPanel = FreqPanelDPOAE('carrier frequency', EXP);

% ---Levels
Levels = SPLpanelDPOAE('-', EXP);
P.DAC='Both';

if ~isequal(EXP.Audio.DAchannelsUsed,P.DAC)
    error('The experiment must be binaural');
end

% ---Durations
Dur = DurPanel('Durations', EXP, '', 'basicsonly');
% ---Pres
Pres = PresentationPanel_XY('F2','L1');
% Pres = PresentationPanel;
% ---Summary
summ = Summary(25);

%====================
Params=GUIpiece('Params'); % upper half of GUI: parameters

Params = add(Params, summ);
Params = add(Params, FPanel, nextto(summ), [10 0]);
Params = add(Params, Levels, nextto(FPanel), [10 0]);
Params = add(Params, Dur, below(FPanel) ,[0 10]);
Params = add(Params, Pres, below(Levels) ,[0 10]);
Params = add(Params, PlayTime, below(Pres) , [-300 10]);




