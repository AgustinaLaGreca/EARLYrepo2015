function iv=isvoid(Q)
% PulldownMenu/isvoid - true for void PulldownMenu object
%   isvoid(P) returns true (1) when P is a void PulldownMenu object, i.e.
%   the result of calling PulldownMenu with no input arguments. If P is an
%   array of PulldownMenu objects, isvoid returns a logical array having
%   the same size.
%
%   See PulldownMenu.


iv = cellfun('isempty', {Q.Name});