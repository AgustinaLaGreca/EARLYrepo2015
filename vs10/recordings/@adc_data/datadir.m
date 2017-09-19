function D = datadir(D);
% adc_data/datadir - full path of local directory containing sampled data
%   datadir(D) returns the full path of the directory where the sampled
%   data of adc_data object D reside. This directory is identified based on
%   two pieces of information: the general data folder returned by 
%   Experiment/parentfolder and the the experiment name as stored in D.
%
%   See also adc_data/samples, experiment/folder.

% this is somewhat of a hack. First reconstruct experiment name
[dum, ExpName] = fileparts(D.Data.Samples.Dir);
% Then abuse Experiment/folder to locate the named experiment
D = folder(dummy(experiment, ExpName));










