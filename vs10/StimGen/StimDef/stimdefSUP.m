function Params = stimdefSup(EXP);
% stimdefSup - definition of stimulus and GUI for SUP stimulus paradigm
%    P=stimdefSUP(EXP) returns the definition for the SUP (suppression)
%    stimulus paradigm. The definition P is a GUIpiece that 
%    can be rendered by GUIpiece/draw. Stimulus definition like stimdefSup 
%    are usually called by StimGUI, which combines the parameter panels with
%    a generic part of stimulus GUIs. The input argument EXP contains 
%    Experiment definition, which co-determines the realization of
%    the stimulus: availability of DAC channels, calibration, recording
%    side, etc.
%
%    See also stimGUI, stimDefDir, Experiment, makestimSup.

% ---Probe sweep
Fsweep = FrequencyStepper('probe frequency', EXP, '', 'notol', 'nobinaural');
%---Supressor tone complex
Sup = local_sup(EXP);
% ---Levels
Levels = SPLpanel('-', EXP, '', 'Probe');
% ---Pres
Pres = PresentationPanel;

%====================
Params = GUIpiece('Params'); % upper half of GUI: parameters
Params = add(Params, Fsweep, 'below', [5 0]);
Params = add(Params, Sup, 'below', [3 10]);
Params = add(Params, Levels, nextto(Fsweep), [30 0]);
Params = add(Params, Pres, below(Levels) ,[85 5]);
Params = add(Params, PlayTime(), below(Pres) , [-55 32]);

%===========================

function T = local_sup(EXP);
if isequal('Both', EXP.AudioChannelsUsed), 
    Nchan = 2;
    PairStr = ' Pairs of numbers are interpreted as [left right].';
else, % single Audio channel
    Nchan = 1;
    PairStr = ''; 
end
% # DAC channels fixes the allowed multiplicity of user-specied numbers
CF = ParamQuery('SupCenterFreq', 'center:', '22987 ', 'Hz', ...
    'rreal/positive', 'Center frequency of suppressor.', 1);
BW = ParamQuery('SupBandWidth', 'width:', '22987 ', 'Hz', ...
    'rreal/positive', 'Bandwidth of suppressor.', 1);
Fbase = ParamQuery('SupBaseFreq', 'base:', '3.75 ', 'Hz', ...
    'rreal/positive', 'Base frequency of suppressor.', 1);
Qlty  = ParamQuery('SupQuality', 'quality:', '', {'perfect' 'approx'}, ...
    'rreal/positive', 'Zwuis quality, i.e., uniqueness of 3rd-order DPs.', 1, 'fontsize', 8);
SPL = ParamQuery('SupSPL', 'SPL:', '-34.6 -45.8', 'dB SPL', ...
    'rreal', ['Intensity of suppressor.' PairStr], Nchan);
Nc = ParamQuery('SupNcomp', '#comp:', '12 ', '', ...
    'posint', 'Number of components of suppressor.', 1);
Dur = ParamQuery('SupDur', 'dur:', '23456 ', 'ms', ...
    'rreal/positive', 'Duration for suppressor+probe. Will be rounded to integer # base cycles.', 1);
OnsetDelay = ParamQuery('OnsetDelay', 'delay:', '1500', 'ms', ...
    'rreal/nonnegative', 'Delay of stimulus onset (common to both DA channels).',1);
Seed = ParamQuery('SupSeed', 'seed:', '234566 ', '', ...
    'rseed', 'Random seed for zwuis realization. Nan = refresh every presentation.', 1);
Msr = messenger('SuppResult', '... overlaps',1, 'ForegroundColor', [0.1 0 0.35], 'fontweight', 'bold');

T = GUIpanel('Suppressor', 'suppressor');
T = add(T, CF, 'below', [0 0]);             T = add(T, Qlty, nextto(CF), [7 4]);
T = add(T, BW, alignedwith(CF), [-1 -5]);   T = add(T, Fbase, nextto(BW), [10 0]);
T = add(T, SPL, alignedwith(BW), [-12 -5]); T = add(T, Nc, nextto(SPL), [10 0]);
T = add(T, Dur, alignedwith(SPL), [0 -5]);  T = add(T, OnsetDelay, nextto(Dur), [10 0]);
T = add(T, Seed, below(Dur), [15 0]);
T = add(T, Msr, nextto(Seed), [10 8]);



