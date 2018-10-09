function C=draw(h, C, XY);
% CycleList/draw - draw cycleList object.
%   draw(h,C) renders cycleList C in the uimenu control with handle h.
%
%   See CycleList, CycleList/refresh, GUIpiece/draw, GUIpiece.

C.parentHandle = h;
% create invisible uimenu that serves as "container"
C.Handle = uimenu(C.parentHandle,...
        'visible', 'off', ...
        'tag', ['Cycle list ' C.Name]);

% Create uimenuitems
Sepa = 'on';
for ii=1:C.Nmax,
    C.ItemHandles(ii) = uimenu(C.parentHandle, ...
        'visible', 'off', ...
        'separator', Sepa, ...
        'Callback', {@select cycleList C.Handle ii}, ... % the first arg serves to invoke the select method
        'Label', ['unused item #' num2str(ii)], ...
        C.uimenuProps);
    Sepa = 'off'; % only first item has separator
end

% Store C in GUIdata of figure
figh = parentfigh(h);
figh = figh.Number;
CL = getGUIdata(figh, 'CycleList', []);
CL = [CL C];
setGUIdata(figh, 'CycleList', CL);

% if no items are present, try to retrieve them from file
if isempty(C.Items), 
    C = load(C);
end

% now synchronize the rendering with the items in C.Items. This is
% delegated to the refresh method, which also updates the GUIdata of the
% figure and the userdata of the "container" iumenu with handle C.Handle.
C=refresh(C);



