function [A, Mess]=criticalDownload(A,St, Verb, Name);
% action/criticalDownload - download Action object and check its status
%    [A,Mess] = criticalDownload(A, 'raw', 'eat') downloads A and test
%    whether Status(A) equals 'cooked'. If A passes the test, Mess is
%    empty. If A fails the test, Mess equals 'Cannot eat a <CCC> object
%    whose status is <XXX>'. Here <CCC> equals class(A) and <XXX> equals
%    status(A). 
%
%    [A,Mess] = criticalDownload(A, {'halal' 'kosher' '5-dimensional'}, 'eat') 
%    allows A to have one of multiple statuses.
%
%    See also Dynamic/download, Action/status.

if nargin<4, Name = ''; end

A = download(A);
[St, N] = cellify(St);
okay = 0;
for ii=1:N,
    okay = okay || isequal(lower(St{ii}), lower(A.Status));
end
if okay,
    Mess = '';
elseif isempty(Name),
    Mess = ['Cannot ' Verb ' a ' class(A) ' object whose status is ''' A.Status  '''.'];
else,
    Mess = ['Cannot ' Verb ' ' class(A) ' object ''' Name '''. Its status is ''' A.Status  '''.'];
end


