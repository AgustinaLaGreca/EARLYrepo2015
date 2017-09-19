function A=prepare(A);
% actabit/prepare - perfrom the prepare action of actabit object
%    A=prepare(A) checks A's todolist for a prepare command; if one is found,
%    it is executed. Thereafter A's status is set to 'prepared'.
%
%    See also actabit.

eval(IamAt);
[A, Mess] = criticalDownload(A, 'initialized', 'prepare');
error(Mess);

if isfield(A.Todolist, mfilename),
    feval(A.Todolist.(mfilename){:});
end

status(A, 'prepared'); % implicit upload here




