function setdefault(P);
% dataviewparam/setdefault - specify the default value for dataviewparam 
%   setdefault(P) sets the default value of type-P dataviewparam objects to
%   P. Dataviewparam belong to the same type if they share the same
%   Dataviewer.
%
%   See also dataviewparam/getdefault, dataviewparam/save.

% P may not be void ...
if isvoid(P),
    error('Cannot set default to void a dataviewparam.');
end
% ... or obsolete
[Template, ParamGUI] = P.Dataviewer(dataset(),'params');
if ~isequal(fieldnames(Template), fieldnames(P.Param)),
    error('An obsolete dataviewparam object cannot serve as default.');
end

CFN = fullfile(GUIdefaultsDir('DataviewParam', 'create'), 'Default.dvparam');
CacheParam = P.Dataviewer;
putcache(CFN, 1e3, CacheParam, P);
save(P, 'def');




