function D=dynamic(Ax,Name);
% dynamic - constructor of dynamic object
%    Dynamic returns an dynamic object having name Foo.
%
%    Dynamic(Ax) returns a dynamic object with access struct Ax as
%    described in GUI access.
%
%    Dynamic(figh,'Foo') is the same as Dynamic(GUIaccess(figh,'Foo')).
%
%    Dynamic is an abstact class, acquiring its usefulness through
%    subclasses like dataset, collectdata, etc.
%    The dynamic class allows objects to be dynamically accessed, as if 
%    they have pointers. To this end, the following methods are defined:
%
%         upload: store yourself
%       download: retrieve yourself
%        offload: remove yourself
%        statify: make read-only
%       hasbeenuploaded: true if previosuly uploaded
%
%    See also GUIaccess, dataset, collectdata.

if nargin<1, Ax = []; end
if nargin<2, Name = ''; end

if nargin>1, Ax = GUIaccess(Ax,Name); end


D.access = Ax;
D.isstatic = 0;
D.name = Name;
D = class(D,mfilename);
if ~isempty(Ax), % upload
    Ax.put(D);
end



