function C=refresh(C);
% CycleList/refresh - refresh rendering of cycleList object.
%   refresh(C) refreshed the rendering of cycleList C. 
%   C must be drawn before it can be refreshed.
%
%   See CycleList, CycleList/draw, GUIpiece/draw, GUIpiece.

% make sure not to attempt to show more items than possible
Nitem = min(C.Nmax, numel(C.Items));
C.Items = C.Items(1:Nitem);
% Update uimenu items of C according to C.Items.
set(C.ItemHandles, 'visible', 'off');
for ii=1:Nitem,
    it = C.Items(ii);
    h = C.ItemHandles(ii);
    set(h, 'Label', ['&' num2str(ii) ' ' it.Label], ...
        'Userdata', it.Userdata, ...
        'visible', 'on');
end

% update userdata of GUI figure
figh = parentfigh(C.Handle);
figh = figh.Number;
CL = getGUIdata(figh, 'CycleList', []);
% replace the correct one in array CL
iam = find(CL, C.Name);
CL(iam) = C;
setGUIdata(figh, 'CycleList', CL);

% update userdata in container uimenu
set(C.Handle, 'userdata', C);

% store current State in file
Mess=save(C);
if ~isempty(Mess),
%     GUImessage(figh, {['Problem saving CycleList ' C.Name], Mess});
    warning(['Problem saving CycleList ' char(10) C.Name char(10) Mess]);
end


