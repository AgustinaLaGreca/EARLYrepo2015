function P=draw(h, P, XY);
% PulldownMenu/draw - draw PulldownMenu object.
%   draw(figh, P) renders PulldownMenu object P in the GUI with handle figh.
%   Items contained in P (see pulldownMenu/additem) are rendered
%   recursively.
%
%   draw(uimenuh, P) renders PulldownMenu object P as an item of an
%   existing uimenu object with handle uimenuh.
%
%   See CycleList, GUIpiece/draw, GUIpiece.

% create topmost uimenu
P.Handle = uimenu(h, 'Label', P.Label, P.uimenuProps);
P.Handle = double(P.Handle);
% draw children ("items") one by one & store their handles
for ii=1:numel(P.Items),
    it = P.Items{ii};
    switch lower(class(it)),
        case  'struct', % elementary item
            ih = uimenu(P.Handle, it); % uimenu with properties in struct it
        case  'pulldownmenu', % another pulldownmenu: recursive
            p = draw(P.Handle, it); 
            ih = p.Handle;
        case  'cyclelist', % another pulldownmenu: recursive
            c = draw(P.Handle, it); 
            ih = c.Handle;
            P.CycleItems = [P.CycleItems c];
    end
    ih = double(ih); %Jan 2018
    P.ItemHandles(ii) = ih;
end

% Store P (including its handles) in userdata of uimenu 
set(P.Handle, 'userdata', P);

% update userdata of GUI figure
figh = parentfigh(P.Handle);
figh = figh.Number;
PD = getGUIdata(figh, 'PulldownMenu', []);
PD = [PD P];
setGUIdata(figh, 'PulldownMenu', PD);



