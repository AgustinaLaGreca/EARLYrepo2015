function [N, T, C]=lastcurrentname(dum); 
% experiment/lastcurrentname - name of most recent experiment
%    lastcurrentname(dumE), where dumE is a dummy experiment arg, returns 
%    name of the most recently Experiment to be checked in. This involves
%    no testing of existence, TimeOut, etc. An empty string '' is returned
%    when no experiment was ever checked in o nthis computer.
%
%    [Name, T, Curr] = lastcurrentname(dumE) also returns the timeout of
%    the experiment (in hours) and the logical Curr which is True when the
%    experiment is still active, False otherwise.
%
%    See also Experiment/checkin, Experiment/exist, Experiment/load.

N = FromSetupFile('experiment', 'currentExperiment', '-default', '');

if nargout>1,
    T = FromSetupFile('experiment', 'TimeOut', '-default', nan);
end
if nargout>2,
    C = false;
    if ~exist(dum, N), return; end
    E = find(experiment, N);
    if isequal('Finished', state(E)), return; end
    if isvoid(E), return; end
    if ~isnan(T),
        Age_hours = (now-datenum(statusmodified(E)))*24; % # hours since last modified
        C = Age_hours<=T; 
    end
end




