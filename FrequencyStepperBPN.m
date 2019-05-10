function Modsweep=FrequencyStepperBPN(T, EXP, Prefix);
% SAMpanel - generic amplitude modulation panel for stimulus GUIs.
%   S=SAMpanel(Title, EXP) returns a GUIpanel M allowing the 
%   user to specify a fixed sinusiodal amplitude modulation to be applied 
%   to all the stimuli of a series. Guipanel S has title Title. EXP is the 
%   experiment definition, from which the number of DAC channels used 
%   (1 or 2) is determined.
%
%   The paramQuery objects contained in F are named: 
%      StartModFreq: starting frequency in Hz
%  StepModFrequency: step in Hz or Octaves (toggle unit)
%   EndModFrequency: end frequency in Hz
%     AdjustModFreq: toggle selecting which of the above params to adjust
%                    in case StepFrequency does not fit exactly.
%    ModFreqTolMode: toggle selecting whether frequencies should be
%                    realized exactly, or whether memory-saving rounding is
%                    allowed.
%           ModDepth: modulation depth in % (100%==full AM)
%      ModStartPhase: starting phase of modulation in Cycles (0=cos)
%
%   SAMpanel is a helper function for stimulus definitions like stimdefFS.
% 
%   M=SAMpanel(Title, EXP, Prefix) prepends the string Prefix
%   to the paramQuery names, e.g. ModFreq -> NoiseModFreq, etc.
%   Default Prefix is '', i.e., no prefix.
%
%   M=SAMpanel(Title, EXP, Prefix, 1) also includes a Theta query. Theta is
%   the modulation angle.
%
%   Use EvalSAMpanel to read the values from the queries.
%
%   See StimGUI, GUIpanel, EvalSAMpanel, stimdefFS.

[Pref, IncTheta, NoTolerance] = arginDefaults('Prefix/IncludeTheta/NoTol', '', 1, 0);
if (nargin < 3), Prefix = Pref; end
if (nargin < 4), IncludeTheta = IncTheta; end
if (nargin < 5), NoTol = NoTolerance; end
    

% # DAC channels fixes the allowed multiplicity of user-specied numbers
if isequal('Both', EXP.AudioChannelsUsed), 
    Nchan = 2;
    PairStr = ' Pairs of numbers are interpreted as [left right].';
else, % single Audio channel
    Nchan = 1;
    PairStr = ''; 
end
ClickStr = ' Click button to select ';

% ---Modulation
Modsweep = GUIpanel('Modsweep', T);
StartModFreq = ParamQuery([Prefix 'StartModFreq'], 'start freq:', '15000.5 15000.5', 'Hz', ...
    'rreal/nonnegative', ['Starting frequency of series.' PairStr], Nchan);
StepModFreq = ParamQuery([Prefix 'StepModFreq'], 'step freq:', '12000', {'Hz' 'Octave'}, ...
    'rreal/positive', ['Frequency step of series.' ClickStr 'step units.'], Nchan);
EndModFreq = ParamQuery('EndModFreq', 'end freq:', '12000.1 12000.1', 'Hz', ...
    'rreal/nonnegative', ['Last frequency of series.' PairStr], Nchan);
AdjustModFreq = ParamQuery([Prefix 'AdjustModFreq'], 'adjust:', '', {'none' 'start' 'step' 'end'}, ...
    '', ['Choose which parameter to adjust when the stepsize does not exactly fit the start & end values.'], 1,'Fontsiz', 8);
Tol = ParamQuery([Prefix 'ModFreqTolMode'], 'acuity:', '', {'economic' 'exact'}, '', [ ...
    'Exact: no rounding applied;', char(10), 'Economic: allow slight (<1 part per 1000), memory-saving rounding of frequencies;'], ...
    1, 'Fontsiz', 8);

% ModDepth = ParamQuery([Prefix 'ModDepth'], ...
%     'depth:', '20.5 33.3', '%', 'rreal/nonnegative', ...
%     ['Modulation depth.' PairStr], Nchan);
% ModStartPhase = ParamQuery([Prefix 'ModStartPhase'], ...
%     'phase:', '-0.25 -0.33', 'Cycle', 'rreal', ...
%     ['Starting phase of modulation. Zero means cosine phase.' PairStr], Nchan);
% Theta = ParamQuery([Prefix 'ModTheta'], ...
%     'theta:', '-0.25 -0.33', 'Cycle', 'rreal', ...
%     ['Modulation angle. Theta = 0 results in AM; Theta=0.25 is QFM; other values give mixed modulation.' PairStr], Nchan);

CutoffSide = ParamQuery('CutoffSide', 'Select which edge to vary:', '', {'L' 'H'}, ...
    '', 'L varies the lower side, making a HPF, and vice versa.');

Modsweep = add(Modsweep, StartModFreq);
Modsweep = add(Modsweep, StepModFreq, alignedwith(StartModFreq));
Modsweep = add(Modsweep, EndModFreq, alignedwith(StepModFreq));
Modsweep = add(Modsweep, AdjustModFreq, nextto(StepModFreq), [10 0]);
Modsweep = add(Modsweep, CutoffSide, below(EndModFreq));
% Modsweep = add(Modsweep,ModDepth, alignedwith(EndModFreq));
% Modsweep = add(Modsweep,ModStartPhase, alignedwith(ModDepth));
% Modsweep = add(Modsweep,ModITD, alignedwith(ModStartPhase));

% if IncludeTheta,
%     Modsweep = add(Modsweep,Theta, 'aligned', [0 -5]);
% end
% if ~isequal('notol', NoTol),
%     Modsweep = add(Modsweep, Tol, alignedwith(AdjustModFreq) , [0 -10]);
end



