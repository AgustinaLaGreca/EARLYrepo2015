function [W, S] = PassYourselfToContextMenu(Src, Ev, VoidQuery);
% ParamQuery/PassYourselfToContextMenu - store ParamQuery in uicontextmenu userdata
%   PassYourselfToContextMenu(Src, Ev, VoidQuery) sets the userdata
%   property of the UIcontextmenu of the calling uicontrol
%   This allows the uimenus of the uicontextmenu to find out which
%   ParamQuery triggered the call.

Q = get(Src, 'userdata');
hct = get(Src,'uicontextmenu');
if isSingleHandle(hct),
    set(hct,'userdata', Q);
end











