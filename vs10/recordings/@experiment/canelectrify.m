function cr = canelectrify(E);
% experiment/canelectrify - True if electrical stim is active for experiment
%   canelectrify(E) - returns True if Experiment object E has any type of
%   electrical stimulation specified, False otherwise. 
%
%   See also Experiment/edit.

cr = ~isequal('-', E.Electra.DA1_target) ...
    || ~isequal('-', E.Electra.DA2_target);







