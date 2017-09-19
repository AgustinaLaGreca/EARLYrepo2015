function ic=iscurrent(E);
% experiment/iscurrent - true if experiment is the current one.
%    iscurrent(E) returns True when E is the current experiment, False
%    otherwise. An experiment is "current" only when at this computer and 
%    at this time the recordings may be added to it. Even when the 
%    experiment is ongoing, it is only current at the measurement computer.
%    Experiments are considered equal when the have the same
%    name (note that it is impossible to make a namesake experiment).
%    By definition a void experiment is not current.
%
%    See also experiment/current, experiment/makecurrent,
%    experiment/finish.

if ~isempty(sys3devicelist) || FromSetupFile('experiment','doHostExperiments', '-default', false),
    [N, T, isC]=lastcurrentname(E);
    ic = isC && isequal(N, name(E));
else,
    ic = false;
end








