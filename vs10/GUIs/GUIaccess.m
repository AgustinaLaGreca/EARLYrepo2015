function S=GUIaccess(figh,FN);
% GUIaccess - access to GUIdata
%    S=GUIacces(figh,'Foo') returns an Address object whose fields are 
%    function handles that get, set, and remove GUIdata, respectively:
%       X=S.get() is equivalent to X = getGUIdata(figh,'Foo')
%       S.put(X) is equivalent to setGUIdata(figh,'Foo',X)
%       S.rm() is equivalent to rmGUIdata(figh,'Foo')
%
%    The names may also be subfields, like 'Foo.fld', etc, as described in
%    getGUIdata.
%
%    See also Address, GUIputter, GUIkiller, GUIacces, getGUIdata.

S = address(GUIgetter(figh,FN), GUIputter(figh,FN), GUIkiller(figh,FN), figh);




