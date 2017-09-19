function wrapup(A);
% action/wrapup - wrapup action of action object
%    A=wrapup(A) wraps up after the action of A is finished
%
%    This generic action/wrapup method only deletes the timer (if any) and 
%    changes A's status to 'wrappedup'. For many subclasses, wrapup will need 
%    overloading. This generic action/wrapup serves as a "catcher" that 
%    prevents a proliferation of timers. It may still be convenient to use
%    end the overloaded XXX/wrapup with the command wrapup(X.action), which
%    calls the generic action/wrapup.
%
%    See also playsound/wrapup, action/clear.

eval(IamAt);
[A, Mess] = criticalDownload(A, {'finished', 'stopped'}, 'wrap up');
error(Mess);

% get rid of timer - don't stop it first , that should have been handled by
% either finishWhenReady or stop.
if ~isempty(A.thandle), 
    delete(A.thandle); 
    A.thandle = [];
end;
% change status & upload
status(A, 'wrappedup');


