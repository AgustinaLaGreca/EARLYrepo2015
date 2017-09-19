function St=status(A,St);
% action/status - get/set status of action object
%    status(A) returns the current status of action object A. If A is 
%    dynamic, it is downloaded before returning its status.
%
%    A=status(A,'Foo') sets the status of A to Foo. If A is dynamic, it is
%    downloaded before, and uploaded after its status has been set.
%
%    See also collectdata, playaudio, createtimer.

HBU = isa(A,'dynamic') && hasbeenuploaded(A);
if HBU, A = download(A); end
if nargin<2, % get
    St = A.Status;
else, % set
    A.Status = St;
    if HBU, upload(A); end
    St = A;
end


