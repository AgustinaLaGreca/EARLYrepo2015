function cr = canrecord(E);
% experiment/canrecord - True if any recording was specified for experiment
%   canrecord(E) - returns True if Experiment object E has any type of
%   recording specified, False otherwise. 
%
%   See also Experiment/edit.

cr = false; % pessimistic default
if ~isvoid(E),
    [RS, FNS] = recordingsources(E);
    for ii=1:numel(FNS),
        if ~isequal('-', RS.(FNS{ii}).DataType),
            cr = true;
            return;
        end
    end
end







