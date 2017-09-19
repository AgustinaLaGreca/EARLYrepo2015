function Params = stimdefMTS(EXP);
% stimdefMTS - definition of stimulus and GUI for MTS stimulus paradigm
%    P=stimdefMTS(EXP) returns the definition for the MTS (multitone suppression)
%    stimulus paradigm. The definition P is a GUIpiece that 
%    can be rendered by GUIpiece/draw. Stimulus definition like stimdefMTS 
%    are usually called by StimGUI, which combines the parameter panels with
%    a generic part of stimulus GUIs. The input argument EXP contains 
%    Experiment definition, which co-determines the realization of
%    the stimulus: availability of DAC channels, calibration, recording
%    side, etc.
%
%    See also stimGUI, stimDefDir, Experiment, makestimMTS.

% ---Frequencies
Tones = local_tonePanel(EXP);
SPLsweep = SPLstepper('varied SPL', EXP);

Probe = FrequencyStepper('probe', EXP, '', 'notol', 'nobinaural');
Lock = ParamQuery('LockProbe', 'lock:', '', {'-' 'to suppr.'}, ...
    '', ...
    {'Click button to have frequencies locked to suppressor frequency or not.' ...
    'Locking to the suppressor means that all probe frequencies are rounded to the' ...
    'nearest multiple of the suppressor frequency.' ...
    'Without locking, the frequencies of probe and suppressor' ...
    'are rounded to multiples of ~5 Hz to increase efficiency.'});
Probe = add(Probe, Lock, 'below lowest', [10 -5]);
%---Supressor
Sup = local_sup(EXP);
%--- Durations
Dur = DurPanel('-', EXP, '', 'basicsonly_mono');
% ---Levels
Levels = SPLpanel('-', EXP, '', 'Probe');
% ---Pres
Pres = PresentationPanel;

%====================
Params = GUIpiece('Params'); % upper half of GUI: parameters
Params = add(Params, Probe, 'below', [5 0]);
Params = add(Params, Dur, below(Probe) ,[0 5]);
Params = add(Params, Levels, nextto(Probe), [30 0]);
Params = add(Params, Sup, below(Levels), [0 0]);
Params = add(Params, Pres, nextto(Sup) ,[15 10]);
Params = add(Params, PlayTime(), below(Dur) , [50 10]);

%===========================
function T = local_tonePanel();
T = GUIpanel('ToneComplex', 'tone complex');
% frequencies typed as an array
Freq = ParamQuery('Freq', 'freq:', '22987 22987 22987 22987 22987 ', 'Hz', ...
    'rreal/positive', {'Whitespace-delimited individual frequencies, max 7 tones' 'The first freq corresponds to the tone whose SPL is varied.'}, 7);
SPL = ParamQuery('fixedSPL', 'SPL:', '-34.6 ', 'dB SPL', ...
    'rreal', 'Intensity of all tones except first one.', 1);
T = add(T, Freq, 'below', [0 0]);
T = add(T, SPL, alignedwith(Freq), [-1 -5]); 



