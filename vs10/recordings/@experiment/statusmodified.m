function n=statusmodified(E);
% experiment/statusmodified - statusmodified data string of Experiment object
%   statusmodified(E) returns the statusmodified property of Experiment 
%   object E.
%
%   See also Experiment/modified

if isvoid(E),
    n = [];
else,
    n = E.Status.StatusModified;
    if isempty(n),
        n = E.ID.Started;
    end
end






