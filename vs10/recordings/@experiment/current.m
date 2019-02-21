function E=current(dumE); 
% experiment/current - current experiment.
%    current(dumE), where dumE is a dummy experiment arg, returns the
%    current experiment. If none is defined, or if the one defined
%    cannot be found, Current returns a void experiment object.
%    The being-current of an experiment is remembered across Matlab
%    sessions, but an TimeOut  (default 6 hours) is employed: if the
%    experiment status has not been modified during his interval, it loses
%    its status as current experiment.
%
%    See also experiment/makecurrent, experiment/touch.

% Marta 11/2018 - Useless line
%if nargin<2, Ndays=1; end % default expiration interval

[Nam, ~, isC] = lastcurrentname(experiment);
if isC && exist(experiment(), Nam)
    E = find(experiment(), Nam);
else
    E = experiment();
end







