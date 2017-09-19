function makecurrent(E, TimeOut);
% experiment/makecurrent - make experiment the current one.
%    makecurrent(E) makes E the current experiment. This role is remembered
%    accross Matlab sessions.
%
%    makecurrent(E, T) uses a TimeOut of T hours. Default TimeOut is 6
%    hours. If the experiment has not been updated for T hours, it is no 
%    longer current.
%
%    See also experiment/iscurrent, experiment/current, experiment/finish.

% checkout current experiment
CE = current(experiment);
if ~isvoid(CE) && ~isequal(name(CE), name(E)),
    finish(CE, ' (makecurrent)');
end

if nargin<2, TimeOut = 6; end; % default 6 hours timeout
ToSetupFile('experiment','-propval','currentExperiment', name(E), 'TimeOut', TimeOut);







