function clear(A);
% action/clear - clear action object
%    clear(A) removes A from history by offloading it.
%    If this is not enough for a certain Action subclass, overload clear.
%
%    See also action/wrapup.

eval(IamAt);
[A, Mess] = criticalDownload(A, 'wrappedup', 'clear');
error(Mess);
offload(A);

