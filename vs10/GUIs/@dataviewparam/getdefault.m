function [P, Mess] = getdefault(dum, Dataviewer);
% dataviewparam/getdefault - get the default value for dataviewparam 
%   P = getdefault(dataviewparam(), @Foo) gets the default 
%   value of dataviewparam objects whose Dataviewer is Foo. If no default 
%   was explicitly defined for the type-Foo dataviewparam objects,
%   P will be void.
%
%   getdefault(P) for non-void P is the same as 
%   getdefault(dataviewparam(), P.Dataviewer)
%
%   See also dataviewparam/setdefault.

% callbacks: void P and extra args

if nargin<2, 
    if isvoid(dum),
        error('Cannot give default for void dataviewparam object w/o specifying the dataviewer.');
    end
    Dataviewer = dum.Dataviewer;
end
Mess = ''; 
if ischar(Dataviewer), Dataviewer = fhandle(Dataviewer); end
% get most recent info on the parameters used by Dataviewer
[Template, ParamGUI] = Dataviewer(dataset(),'params');
ParamNames = fieldnames(Template);
CFN = fullfile(GUIdefaultsDir, 'DataviewParam', 'Default.dvparam');
CacheParam = fhandle(Dataviewer); % allow fcn names, but convert them
P = getcache(CFN, CacheParam);
needInitialization = 0;
if isempty(P), % no default was ever defined or it has been erased.
    P = dum;
    P.Dataviewer = Dataviewer;
    needInitialization = 1;
elseif  ~isequal(ParamNames, fieldnames(P.Param)), % parameter collection has changed since last storage of default
    needInitialization = 1;
else, % merge with newest info from Dataviewer
    P.Param = union(Template, P.Param);
    P.Param = structpart(P.Param, fieldnames(Template));
end

if needInitialization,
    P = local_init(P, Template, ParamGUI);
    if isvoid(P), 
        Mess = 'Unable to create dataviewparam: default values are lacking or obsolete.';
    end
end
if nargout<2, error(Mess); end
%===============
function P = local_init(P, Template, GUI);
uiwait(msgbox({['Default display/analysis parameters for dataset/' char(P.Dataviewer) ' '], ...
    'have not been specified or are obsolete.', 'Please provide them using the following GUI'}, ...
    ['Initialization of ' char(P.Dataviewer) ' parameters.'], 'modal'))
P = edit(P, 'Provide default values.');
if ~isvoid(P), % i.e., if user didn't cancel
    setdefault(P); % save as default
end




