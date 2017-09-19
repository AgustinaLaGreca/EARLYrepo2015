function G=GUIgetter(figh,FN);
% GUIgetter - pointer to get GUIdata
%    PTR=GUIgetter(figh,'Foo') returns a function handle ("pointer") PTR
%    that provides passive access to the 'Foo' field of GUIdata. That is, 
%    X = PTR() is equivalent to getGUIdata(figh,'Foo'). 
%
%    See also GUIputter, GUIkiller, GUIacces, getGUIdata.

G = @()getGUIdata(figh,FN);




