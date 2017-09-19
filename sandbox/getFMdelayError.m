function [ Delay ] = getFMdelayError( experiment_name, dataset_id)
% Input arguments:
%   experiment_name: Name of the experiment
%   dataset_id: id number of the dataset can be found in the first column
%               of databrowse fot this experiment
%
% Output arguments:
%   Delay: Error in delay in miliseconds
%
% Created by Adriaan Lambrechts

ds=read(dataset,experiment_name,dataset_id);

[~, ~, ~, Delay]=calibrate(ds.ID.Experiment, ds.Stim.Fsam, ds.Stim.DAC(1),[]);

end

