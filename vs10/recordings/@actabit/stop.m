function A=stop(A);
% actabit/stop - perform the stop action of an actabit object
%    A=stop(A) checks A's todo list for a stop command; if one is found,
%    it is executed. Thereafter A's status is set to 'stopped'.
%    Note that A, once started, has status 'finished'. Therefore, actabit
%    requires A to have status 'finished' on entry.
%
%    To understand why it is necessary to change the Status from 'finished'
%    to 'stopped', see the comments in the local_allSilent function in
%    GUIaction.
%
%    See also actabit, GUIaction.

eval(IamAt);
%whoiscalling
[A, Mess] = criticalDownload(A, 'finished', 'stop');
if ~isempty(Mess),
    A
end
error(Mess);

if isfield(A.Todolist, mfilename),
    feval(A.Todolist.(mfilename){:});
end

status(A, 'stopped'); % implicit upload here



