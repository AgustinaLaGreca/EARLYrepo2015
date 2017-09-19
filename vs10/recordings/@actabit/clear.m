function A=clear(A);
% actabit/clear - perfrom the clear action of actabit object
%    A=clear(A) checks A's todolist for a clear command; if one is found,
%    it is executed. Thereafter A is offloaded.
%
%    See also actabit.

eval(IamAt);
[A, Mess] = criticalDownload(A, 'wrappedup', 'clear');
error(Mess);

if isfield(A.Todolist, mfilename),
    feval(A.Todolist.(mfilename){:});
end

offload(A);



