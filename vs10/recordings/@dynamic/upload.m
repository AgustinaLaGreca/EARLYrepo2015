function [D,Ax] = upload(D, Ax, Name);
% Dynamic/upload - upload current version of dynamic object
%
%   D=upload(D, Ax) stores the current value of D in the address Ax. 
%   Ax is a struct with fields get, put and rm as described in GUIaccess.
%
%   [D,Ax]=upload(D, figh, 'Foo') stores D in the GUI with figure handle figh
%   under the name Foo. An Ax struct (see above) is returned.
%
%   upload(D) assumes that D has been uploaded before and uploads D to its
%   previously defined address.
%
%   In all uses of upload, if any previosuly uploaded versions of D are 
%   removed, even if they were stored in a different address.
%
%   Updated version of D can be retrieved by download(D) performed on a
%   previous version of D. This is a hack to use pointer-like programming
%   constructs. 
%   
%   See Dynamic/download, Dynamic/offload, GUIaccess.

if nargin<2, Ax = D.access; end
if nargin<3, Name = ''; end

if isempty(Ax),
    error('Missing or void upload destination of dynamic object.');
elseif D.isstatic && ~isequal('grace', Name),
    error('Attempt to upload read-only (statified) dynamic object.');
elseif ishandle(Ax), % GUI handle
    Ax = GUIaccess(Ax, Name);
end
D.access = Ax;
if ~isempty(Name), D.name = Name; end
put(Ax,D);


