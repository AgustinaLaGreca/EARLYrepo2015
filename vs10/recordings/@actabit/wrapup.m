function A=wrapup(A);
% actabit/wrapup - perform the wrapup action of an actabit object
%    A=start(A) checks A's todolist for wrapup and wrapupstop calls.
%    The action depends on A's status. If status(A) is 'finished', wrapup
%    invokes the 'wrapup' call (if present). If status(A) is 'stopped',
%    the 'wrapupstopped' call is invoked when it exists. If it does not
%    exist, the 'wrapup' call is invoked if that exists.
%
%    See also actabit.

eval(IamAt);
[A, Mess] = criticalDownload(A, {'finished' 'stopped'}, 'wrapup');
error(Mess);

% if status is stopped, wrapupstopped has priority if it exists
if isequal('stopped', status(A)) && isfield(A.Todolist, 'wrapupstopped'),
    feval(A.Todolist.wrapupstopped{:});
elseif isfield(A.Todolist, 'wrapup'), % regular 'wrapup' call is default
    feval(A.Todolist.wrapup{:});
end

status(A, 'wrappedup'); % implicit upload here



