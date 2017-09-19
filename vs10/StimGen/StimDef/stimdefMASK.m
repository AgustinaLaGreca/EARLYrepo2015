function Params = stimdefMASK(EXP);
% stimdefMASK - definition of stimulus and GUI for MASK stimulus paradigm
%    P=stimdefMASK(EXP) returns the definition for the MASK (tone in noise)
%    stimulus paradigm. The definition P is a GUIpiece that 
%    can be rendered by GUIpiece/draw. Stimulus definition like stimdefNPHI 
%    are usually called by StimGUI, which combines the parameter panels with
%    a generic part of stimulus GUIs. The input argument EXP contains 
%    Experiment definition, which co-determines the realization of
%    the stimulus: availability of DAC channels, calibration, recording
%    side, etc.
%
%    See also stimGUI, stimDefDir, Experiment, makestimMASK.

% ---Warp sweep
WarpSweep = local_warp;
%---Masked tone
Tone = local_tone(EXP);
% ---Durations
Dur = DurPanel('Noise timing', EXP, 'Mnoise'); 
% ---Pres
Pres = PresentationPanel;
% ---Noise
Noise = NoisePanel('noise param', EXP, 'Mnoise');
%====================
Params = GUIpiece('Params'); % upper half of GUI: parameters
Params = add(Params, WarpSweep, 'below', [5 0]);
Params = add(Params, Tone, 'below', [12 10]);
Params = add(Params, Noise, nextto(WarpSweep), [10 0]);
Params = add(Params, Dur, below(Noise) ,[0 3]);
Params = add(Params, Pres, nextto(Dur) ,[5 -2]);
Params = add(Params, PlayTime(), below(Dur) , [0 12]);

%===========================
function W = local_warp;
StartWarp = ParamQuery('StartWarp', 'start:', '-1.66 ', 'octave', ...
    'rreal', 'Start value of spectral warp.',1);
StepWarp = ParamQuery('StepWarp', 'step:', '0.0001', 'octave', ...
    'rreal/positive', 'Warp step size.', 1);
EndWarp = ParamQuery('EndWarp', 'end:', '-0.5679 ', 'octave', ...
    'rreal', 'End value of spectral warp.', 1);
AdjustWarp = ParamQuery('AdjustWarp', 'adjust:', '', {'none' 'start' 'step' 'end'}, ...
    '', 'Choose which parameter to adjust when the stepsize does not exactly fit the start & end values.', 1,'Fontsiz', 8);
W = GUIpanel('WarpSweep', 'warp sweep');
W = add(W, StartWarp);
W = add(W, StepWarp, alignedwith(StartWarp), [0 -5]);
W = add(W, EndWarp, alignedwith(StepWarp), [0 -5]);
W = add(W, AdjustWarp, nextto(StepWarp), [10 -5]);

function T = local_tone(EXP);
% # DAC channels fixes the allowed multiplicity of user-specied numbers
if isvoid(EXP) || isequal('Both', EXP.AudioChannelsUsed), 
    Nchan = 2;
    PairStr = ' Pairs of numbers are interpreted as [left right].';
else, % single Audio channel
    Nchan = 1;
    PairStr = ''; 
end
Freq = ParamQuery('ToneFreq', 'freq:', '22987 87662 ', 'Hz', ...
    'rreal/positive', ['Unwarped frequency of tone.' PairStr], Nchan);
SPL = ParamQuery('ToneSPL', 'S/N0:', '-34.6  -34.6', 'dB', ...
    'rreal', ['Intensity of tone re spectrum level of noise.' PairStr], Nchan);
ToneOnset = ParamQuery('ToneOnset', 'onset:', '-3456 -3456 ', 'ms', ...
    'rreal', ['Onset of tone relative to noise onset. Negative value -> tone leads noise.' PairStr], Nchan);
ToneDur = ParamQuery('ToneDur', 'dur:', '23456.2 23456.2 ', 'ms', ...
    'rreal/positive', ['Duration of tone including ramps.' PairStr], Nchan);
ToneRise = ParamQuery('ToneRiseDur', 'rise:', '3.2 3.2 ', '', ...
    'rreal/nonnegative', ['Duration of rising ramp.' PairStr], Nchan);
ToneFall = ParamQuery('ToneFallDur', 'fall:', '3.2 3.3', 'ms', ...
    'rreal/nonnegative', ['Duration of falling ramp.' PairStr], Nchan);
T = GUIpanel('Tone', 'tone');
T = add(T, Freq, 'below', [10 0]);
T = add(T, SPL, alignedwith(Freq), [0 -5]);
T = add(T, ToneOnset, alignedwith(SPL), [0 -5]);
T = add(T, ToneDur, alignedwith(ToneOnset), [0 -5]);
T = add(T, ToneRise, alignedwith(ToneDur), [-8 -5]);
T = add(T, ToneFall, nextto(ToneRise), [-10 0]);




