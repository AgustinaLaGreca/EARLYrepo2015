function A=start(A);
% actabit/start - perform the start action of actabit object
%    A=start(A) checks A's todolist for a start command; if one is found,
%    it is executed. Thereafter A's status is set to 'finished'.
%
%    See also actabit.

eval(IamAt);
[A, Mess] = criticalDownload(A, 'prepared', 'start');
error(Mess);

if isfield(A.Todolist, mfilename),
    feval(A.Todolist.(mfilename){:});
end

status(A, 'finished'); % implicit upload here



