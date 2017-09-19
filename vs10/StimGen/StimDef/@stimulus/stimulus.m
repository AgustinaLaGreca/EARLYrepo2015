function S = stimulus(StimType, GUIcolor, Experiment);
% stimulus - construct Stimulus object
%   Stimulus(StimName, GUIcolor) is a Stimulus object. Stimulus is an 
%   abstract class. Its concrete subclasses should have the following 
%   methods
%        define: define GUI and stimulus parameters
%      generate: compute waveforms
%
%   See also stimFS.


[StimType, GUIcolor, Experiment] = arginDefaults('StimType/GUIcolor/Experiment', []);

[Grab, Handle] = deal([]);
GUI = CollectInStruct(GUIcolor, Grab, Handle);
[Param, Fsam, SecParam, Presentation, Waveforms, Attenuation] = deal([]);


S = CollectInStruct(StimType, GUI, Param, Experiment, Fsam, SecParam, ...
    Presentation, Waveforms, Attenuation);
S = class(S, mfilename);







