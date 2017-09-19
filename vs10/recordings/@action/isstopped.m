function is=isstopped(A);
% action/isstopped - true for stopped action
%    isstopped(A) returns True when A's status is 'stopped', False
%    otherwise.
%
%    See also action/isready, action/stop, action/isready.

A = download(A, 'sloppy');
is = isequal('stopped', A.Status);




