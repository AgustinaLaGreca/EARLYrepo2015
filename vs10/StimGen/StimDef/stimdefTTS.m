function Params = stimdefTTS(EXP);
% stimdefTTS - definition of stimulus and GUI for TTS stimulus paradigm
%    P=stimdefTTS(EXP) returns the definition for the TTS (2 tone suppression)
%    stimulus paradigm. The definition P is a GUIpiece that 
%    can be rendered by GUIpiece/draw. Stimulus definition like stimdefTTS 
%    are usually called by StimGUI, which combines the parameter panels with
%    a generic part of stimulus GUIs. The input argument EXP contains 
%    Experiment definition, which co-determines the realization of
%    the stimulus: availability of DAC channels, calibration, recording
%    side, etc.
%
%    See also stimGUI, stimDefDir, Experiment, makestimTTS.

% ---Probe 
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

function T = local_sup(EXP);
if isequal('Both', EXP.AudioChannelsUsed), 
    Nchan = 2;
    PairStr = ' Pairs of numbers are interpreted as [left right].';
else, % single Audio channel
    Nchan = 1;
    PairStr = ''; 
end

T = GUIpanel('Suppressor', 'suppressor');
% # DAC channels fixes the allowed multiplicity of user-specied numbers
SupFreq = ParamQuery('SupFreq', 'freq:', '22987 ', 'Hz', ...
    'rreal/positive', 'Frequency of suppressor.', 1);
SPL = ParamQuery('SupSPL', 'SPL:', '-34.6 -45.8', 'dB SPL', ...
    'rreal', ['Intensity of suppressor.' PairStr], Nchan);
T = add(T, SupFreq, 'below', [0 0]);
T = add(T, SPL, alignedwith(SupFreq), [-1 -5]); 



