function P=ZUREKecho(T, EXP, Prefix)%, haveTypes);
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

%==========ITD GUIpanel=====================
P = GUIpanel('echo', T);
% ITDpref = preferences(EXP, 'ITDconvention');
echoITD = ParamQuery('echoITDfactor', 'ITD factor:', '-10.170', '', 'rreal', ['-1 (or) 1 means opposite of (or) same as the source respectively' ...
                        newline '0 means no ITD' newline 'Other numbers mean increased/decreased by a factor of that number with respect to the source'...
                        newline 'Should be within [-2 2]'], 1);
echoILD = ParamQuery('echoILDfactor', 'ILD factor:', '-10.170', '', 'rreal', ['ILD for the echo' newline ...
    'Same convention as ITD' newline 'Should be within [-2 2]' ], 1);                 
echoSPL = ParamQuery('echoSPL', 'SPL:', '0.170', 'dB', 'rreal/positive', 'SPL of the echo', 1);
echoSeed = ParamQuery([Prefix 'echoSeed'], 'Noise seed:', '844596300', '', ...
    'rseed', 'Random seed used for realization of echo waveform.',1);


Delay = ParamQuery('Delay', 'Lead-Lag Delay :', '0.170', 'ms', 'rreal', 'Lead - Lag delay', 1);

P = add(P, echoITD, 'below', [50 0]);
P = add(P, echoILD, alignedwith(echoITD));
P = add(P, echoSPL, alignedwith(echoILD));
P = add(P, echoSeed, alignedwith(echoSPL));
P = add(P, Delay, 'aligned');
end



