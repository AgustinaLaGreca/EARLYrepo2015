function p=preferences(E, FN);
% experiment/preferences - preferences of Experiment object
%   preferences(E) returns Preferences struct of experiment E.
%
%   preferences(E, 'foo') returns preference named foo.

p = E.Preferences;
if nargin>1,
    p = p.(FN);
end







