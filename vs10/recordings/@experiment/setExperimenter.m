function [ Exp ] = setExperimenter( Exp, newExperimentName )
%Change or Set the name of the Experiment Exp to newExperimentName
%   Detailed explanation goes here
if nargin == 1 && ischar(newExperimentName)
   Exp.ID.Experimenter = newExperimentName;
end


end

