function n=name(E);
% experiment/name - name of Experiment object
%   Name(E) returns the name of Experiment object E.

n = E.ID.Name;
if isempty(n),
    n = '';
end






