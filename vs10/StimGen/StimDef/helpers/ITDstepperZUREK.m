function P=ITDstepper(T, EXP, Prefix)%, haveTypes);
% ITDstepper - generic ITD stepper panel for stimulus GUIs.
%   F=ITDstepper(Title, EXP) returns a GUIpanel F allowing the 
%   user to specify a series of ITDs, using a linear spacing.  
%   The Guipanel F has title Title. EXP is the experiment 
%   definition, from which the number of DAC channels used (1 or 2) is
%   determined.
%
%   The paramQuery objects contained in F are named: 
%         StartITD: starting frequency in Hz
%     StepITD: step in Hz or Octaves (toggle unit)
%      EndITD: end frequency in Hz
%        AdjustITD: toggle selecting which of the above params to adjust
%                    in case StepFrequency does not fit exactly.
%   ITDStepper is a helper function for stimulus definers like 
%   stimdefBINZW.
% 
%   F=ITDstepper(Title, Exp, Prefix) prepends the string Prefix
%   to the paramQuery names, e.g. StartFreq -> ModStartFreq, etc.
%
%  ITDstepper(Title, Exp, Prefix, haveTypes) specifies whether to allow
%  different types of ITD (ongoing, total waveform, etc). 
%  The default is haveTypes = false.
%
%   Use EvalITDstepper to read the values from the queries and to
%   compute the actual frequencies specified by the above step parameters.
%
%   See StimGUI, GUIpanel, EvalFrequencyStepper, makestimFS.

[Prefix, haveTypes] = arginDefaults('Prefix, haveTypes', '', false);

ClickStr = ' Click button to select ';
ITDstr = ['(positive = ' upper(strrep(EXP.ITDconvention, 'Lead', ' Lead')) ')'];

%Levels = GUIpanel('Levels', T);
DACstr = getDACstr(EXP.AudioChannelsUsed, EXP.Recordingside);

%==========ITD GUIpanel=====================
P = GUIpanel('ITD', T);
% ITDpref = preferences(EXP, 'ITDconvention');
startITD = ParamQuery('startITD', 'ITD start:', '-10.170', 'ms', 'rreal', ['Start value of ITD ' ITDstr], 1);
endITD = ParamQuery('endITD', 'ITD end:', '-10.170', 'ms', 'rreal', ['End value of ITD ' ITDstr], 1);
stepITD = ParamQuery('stepITD', 'ITD step:', '0.170', 'ms', 'rreal/positive', 'Step value of ITD ', 1);
ILD = ParamQuery('ILD', 'ILD:', '-10.170', 'dB', 'rreal', ['ILD for the source' newline ...
    '+x increases left ear level by x/2 and decreases right ear level by x/2' newline '-x decreases left ear level by x/2 and increases right ear level by x/2' ], 1);
adjustITD = ParamQuery([Prefix 'AdjustITD'], 'adjust:', '', {'none' 'start' 'step' 'end'}, ...
    '', ['Choose which parameter to adjust when the stepsize does not exactly fit the start & end values.'], 1,'Fontsiz', 8);
% ITDtype = ParamQuery([Prefix 'ITDtype'], 'impose on', '', {'waveform' 'fine' 'gate' 'mod' 'fine+gate' 'fine+mode' 'gate+mod'}, '', ...
%     ['Implementation of ITD. Click to toggle between options.' char(10) ...
%     '    waveform = whole waveform delay' char(10) ...
%     '      gating = delayed gating imposed on nondelayed waveform' char(10)  ...
%     '     ongoing = nondelayed gating imposed on delayed waveform.' ]);

LowFreq = ParamQuery([Prefix 'LowFreq'], ...
    'low:', '1100.1 1100.1', 'Hz', 'rreal/nonnegative', ...
    ['Low cutoff frequency.' ]);
HighFreq = ParamQuery([Prefix 'HighFreq'], ...
   'high:', '1100.1 1100.1', 'Hz', 'rreal/nonnegative', ...
    ['High cutoff frequency.' ]);
Corr = ParamQuery([Prefix 'Corr'], 'corr:', '-0.9997 ', {'I', 'C'}, ...
    'rreal', ['Interaural noise correlation (number between -1 and 1)', char(10), ... 
    'Click button to change the "varied channel" (where mixing is done).'],1);
NoiseSeed = ParamQuery([Prefix 'ConstNoiseSeed'], 'seed:', '844596300', '', ...
    'rseed', 'Random seed used for realization of noise waveform. Specify NaN to refresh seed upon each realization.',1);
SPL = ParamQuery([Prefix 'SPL'], 'level:', '120.5 120.5', {'dB SPL' 'dB/Hz'}, ...
    'rreal', ['Intensity. Click button to switch between overall level (dB SPL) and spectrum level (dB/Hz).' ]);
AltLevel = messenger([Prefix 'AltLevel'], 'S=XXX dB/Hz', 1);
MaxSPL=messenger([Prefix 'MaxSPL'], 'max [**** ****] dB SPL @ [***** *****] Hz    ',1);
    
DAC = ParamQuery('DAC', 'DAC:', '', DACstr, ...
      '', ['Active D/A channels.' ClickStr 'channel(s).']);

P = add(P, startITD, 'below', [0 0]);
P = add(P, stepITD, alignedwith(startITD));
P = add(P, endITD, alignedwith(stepITD));
P = add(P, adjustITD, nextto(stepITD));
P = add(P, ILD, alignedwith(endITD));

% For the noise segment
P = add(P,LowFreq, alignedwith(ILD), [0 0]);
P = add(P,HighFreq, 'aligned', [0 -7]);
P = add(P,SPL, alignedwith(HighFreq), [0 0]);
P = add(P,Corr, 'aligned', [0 -7]);
P = add(P,NoiseSeed, nextto(Corr), [0 0]);
P = add(P,AltLevel, nextto(SPL), [4 10]);
P = add(P, DAC, alignedwith(Corr));
P = add(P,MaxSPL,below(DAC),[17 0]);

% if haveTypes,
%     P = add(P, ITDtype, below(endITD));
end



% P = GUIpanel('Fsweep', T);
% StartFreq = paramquery([Prefix 'StartFreq'], 'start:', '15000.5 15000.5', 'Hz', ...
%     'rreal/positive', ['Starting frequency of series.' PairStr], Nchan);
% StepFreq = paramquery([Prefix 'StepFreq'], 'step:', '12000', {'Hz' 'Octave'}, ...
%     'rreal/positive', ['Frequency step of series.' ClickStr 'step units.'], Nchan);
% EndFreq = paramquery('EndFreq', 'end:', '12000.1 12000.1', 'Hz', ...
%     'rreal/positive', ['Last frequency of series.' PairStr], Nchan);
% AdjustFreq = paramquery([Prefix 'AdjustFreq'], 'adjust:', '', {'none' 'start' 'step' 'end'}, ...
%     '', ['Choose which parameter to adjust when the stepsize does not exactly fit the start & end values.'], 1,'Fontsiz', 8);
% Tol = paramquery([Prefix 'FreqTolMode'], 'acuity:', '', {'economic' 'exact'}, '', [ ...
%     'Exact: no rounding applied;', char(10), 'Economic: allow slight (<1 part per 1000), memory-saving rounding of frequencies;'], ...
%     1, 'Fontsiz', 8);
% Fsweep = add(Fsweep, StartFreq);
% Fsweep = add(Fsweep, StepFreq, alignedwith(StartFreq));
% Fsweep = add(Fsweep, EndFreq, alignedwith(StepFreq));
% Fsweep = add(Fsweep, AdjustFreq, nextto(StepFreq), [10 0]);
% if ~isequal('notol', Flag),
%     Fsweep = add(Fsweep, Tol, alignedwith(AdjustFreq) , [0 -10]);
% end




