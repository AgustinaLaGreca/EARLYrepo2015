function ST = stimFS(Experiment);
% stimFS - construct stimFS object
%    stimFS(Exp) constructs a stimFS object. StimFS is a subclass of Stimulus.
%    Exp is the experiment in which the stimulus is used.
%    The true definition of the stimulus are provided by the following
%    methods
%       define: 
%     generate:

Experiment = arginDefaults('Experiment',[]);

GUIcolor = [0.65 0.69 0.68];
StimType = 'fs';

ST = struct();
ST = class(ST, mfilename, stimulus(StimType, GUIcolor, Experiment));


