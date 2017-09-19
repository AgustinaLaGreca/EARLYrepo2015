function CreateQueryContextMenu(figh,dum);
% ParamQuery/CreateQueryContextMenu - create standard uicontextmenu for ParamQueries
%   CreateQueryContextMenu(figh, Qdum) creates a standard uicontextmenu in 
%   figure figh, which serves as a default context menu for paramqueries. 
%   The uicontextmenu is retrievable by its tag 'QueryContextMenu'.
%
%   The dummy argument Qdum is a void ParamQuery that allows
%   CreateQueryContextMenu to be a method rather than a function.
%
%   See also PassYourselfToContextMenu.

uc=uicontextmenu('parent', figh, 'tag', 'QueryContextMenu');
uimenu(uc, 'Label', 'freeze', 'callback', {@freeze ParamQuery 1});
uimenu(uc, 'Label', 'thaw', 'callback', {@freeze ParamQuery 0});






