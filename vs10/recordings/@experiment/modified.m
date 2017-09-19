function n = modified(E);
% experiment/modified - modified data string of Experiment object
%   modified(E) returns the 'modified' property of Experiment object
%   E, which indicates when E's definition was last changed.
%
%   See also Experiment/statusmodified

if isvoid(E),
    n = [];
    return;
end

n = E.Status.Modified;
if isempty(n),
    n = E.ID.Started;
end






