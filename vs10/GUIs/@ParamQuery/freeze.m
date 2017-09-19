function [W, S] = Freeze(Src, Ev, VoidQuery, DoFreeze);
% ParamQuery/Freeze - effectively disable uicontrol of rendered ParamQuery
%   Freeze(Src, Ev, VoidQuery) freezes the edit and/or pushbutton controls
%   associated with a rendered ParamQuery Q. Q is not directly passed to
%   Freeze, since it is called via a context menu. The convention is that
%   the calling menu item, whose graphics handle is Src, has as parent
%   a context menu whos handle is, say, Uh, who has the rendered
%   ParamQuery as its userdata. VoidQuery is a void ParamQuery object
%   only passed in order to reach the method Freeze.
%
%   Freeze(Src, Ev, VoidQuery,1) is the same as Freeze(Src, Ev, VoidQuery).
%   Freeze(Src, Ev, VoidQuery,0) unfreezes the controls,

% retrieve the rendered paramquery & the handles of its rendering
uh = get(Src, 'parent'); % uicontext menu handle
Q = get(uh,'userdata');
Qh = Q.uiHandles;

if DoFreeze,
    set(Qh.Edit,'enable', 'off');
    if isfield(Qh,'Unit'),
        set(Qh.Unit, 'enable', 'off');
    end
else,
    set(Qh.Edit,'enable', 'on');
    if isfield(Qh,'Unit'),
        set(Qh.Unit, 'enable', 'on');
    end
end














