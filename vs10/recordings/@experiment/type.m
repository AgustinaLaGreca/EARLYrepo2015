function n=type(E);
% experiment/type - Type of experiment object
%   Type(E) returns the type of Experiment object E.
n = E.ID.Type;
if isempty(n), n = ''; end







