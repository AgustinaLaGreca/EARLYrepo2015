function iv=isvoid(E);
% experiment/isvoid - True for void Experiment object
%   isvoid(E) - returns True for if Experiment object E is void, False
%   otherwise.
%
%   See also Experiment/checkin.

iv = isempty(name(E)) || isempty(experimenter(E));







