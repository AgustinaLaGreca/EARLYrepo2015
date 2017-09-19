function A=prepare(A);
% action/prepare - prepare action object for the action
%    A=prepare(A) prepares A for the action it does for a living
%
%    The generic action/prepare only chekcs whether A' status is 
%    'initialized' and then changes it to 'prepared'. Overload if real
%    preparation is needed.
%
%    See also playsound/prepare.

eval(IamAt);
[A, Mess] = criticalDownload(A, 'initialized', 'prepare');
error(Mess);
A = status(A, 'prepared'); % note: uploading is implicit in this call









