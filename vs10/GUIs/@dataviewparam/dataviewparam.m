function [P, Mess]=dataviewparam(Dataviewer, FN);
% dataviewparam - construct dataviewparam object
%  P=Dataviewparam(@Foo) returns the default dataview parameter object for
%  dataset method Foo (the "dataviewer"). In order for this to work, the 
%  dataviewer Foo must return a template for the parameters in struct T and
%  the GUI G for specifying the parameters when called using the syntax 
%     [T, G] = Foo(dataset(), 'params').
%  If the dataviewparam of a given dataviewer has never been used yet, the
%  GUI will be launched when first calling dataviewparam, and the user must 
%  provide the default values.
%
%  See dataset/dotraster for an example.
%
%  [P, Mess] = Dataviewparam('Foo', FN) attempts to read the parameters 
%  from file with name FN. If FN is not a full filename, its location is
%  determined by the dataviewparam/filename method. If the requested 
%  dataviewparam object cannot be retrieved from file, a void dataviewparam 
%  object is returned and Mess contains an appropriate error message.
%
%  The parameters of a dataviewparam object P may be accessed using the 
%  dataviewparam/get method. P can be passed as a 3rd argument to Foo, and
%  edited by the Edit method.
%
%  Dataviewparam() returns a void dataviewparam object.
%
%  See also dotraster, dataset, "methodshelp dataviewparam".

if nargin<2, FN = 'default'; end 
if nargin<1, 
    Dataviewer = @nope; 
    GUIgrab = [];
    FN = '..void..';
else, % % get the most recent versions the parameter collection & its GUI directly from the Dataviewer
    Dataviewer = fhandle(Dataviewer);
    [Template, ParamGUI] = Dataviewer(dataset(), 'params');
end 

Param = []; Mess = '';

if isempty(FN),
    FN = 'default'; % return default values of parameters
end

if ~ischar(FN),
    error('Second input argument must be char string.');
end

if ~isequal(@nope, Dataviewer), % check if Dataviewer exist
    DVname = func2str(Dataviewer);
    if ~ismethod(dataset(), DVname),
        error(['First argument ''' DVname ''' is not a Dataset method.']);
    end
end

switch lower(FN),
    case '..void..',
    case 'default',
        [P, Mess] = getdefault(dataviewparam(), Dataviewer);
        error(Mess);
        return;
    otherwise, % filename
end

P = CollectInStruct(Dataviewer, Param, GUIgrab);
P = class(P, mfilename);





