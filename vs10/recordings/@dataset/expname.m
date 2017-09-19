function N = expname(D)
% Dataset/expname - experiment name of dataset object.
%   Expname(D) or D.ExpName returns the name of the experiment to which D
%   belongs.
%
%   See Dataset, Experiment.

if isa(D, 'pooled_dataset'),
    D = members(D);
    D = D(1);
end
if isvoid(D),
    error('Attempt to obtain experiment name from void dataset.');
end
N = name(D.ID.Experiment);











