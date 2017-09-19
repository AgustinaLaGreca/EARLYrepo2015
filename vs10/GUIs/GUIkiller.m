function P=GUIkiller(figh,FN);
% GUIkiller - pointer to delete GUIdata
%    PTR=GUIkiller(figh,'Foo') returns a function handle ("pointer") PTR
%    that provides fatal access to the 'Foo' field of GUIdata. That is, 
%    PTR() is equivalent to rmGUIdata(figh,'Foo'). 
%
%    See also GUIgetter, GUIkiller, getGUIdata.

P = @()rmGUIdata(figh,FN);




