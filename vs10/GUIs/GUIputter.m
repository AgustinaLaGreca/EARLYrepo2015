function P=GUIputter(figh,FN);
% GUIputter - pointer to set GUIdata
%    PTR=GUIgetter(figh,'Foo') returns a function handle ("pointer") PTR
%    that provides active access to the 'Foo' field of GUIdata. That is, 
%    PTR(X) is equivalent to setGUIdata(figh,'Foo',X). 
%
%    See also GUIgetter, GUIkiller, getGUIdata.

P = @(X)setGUIdata(figh,FN,X);




