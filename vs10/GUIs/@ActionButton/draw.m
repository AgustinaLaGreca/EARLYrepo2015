function Q=draw(hp, Q, XY);
% ActionButton/draw - draw ActionButton object.
%   draw(hp, B, X,Y) draws ActionButton object B in uipanel with handle hp at
%   position (X,Y) in pixels. This function is typically recursively called 
%   by GUIpiece/draw.
%
%   See ActionButton, GUIpiece/draw, GUIpiece.

cpos=get(hp,'position');
X=XY(1); Y=cpos(4)-XY(2)-Q.Extent(2);

% callback in Q is either function handle or char string. Which it is affects
% the passing of additional Left/Right clickside arg.
if ischar(Q.Callback),
    CBleft = Q.Callback;
    CBrite = Q.Callback;
elseif isfhandle(Q.Callback);
    CBleft = {Q.Callback 'Left'};
    CBrite = {Q.Callback 'Right'};
else, error(['ActionButton callback is a ' class(Q.Callback) ', but should be either char string or function handle. See help text of ActionButton.m']);
end
% In case of a toggling ActionButton, wrap the callback fcn in a call to a
% standard callback.
if iscellstr(Q.String),
    CBleft = {@actionToggleCB ActionButton CBleft}; % first arg is dummy ActionButton object needed to ...
    CBrite = {@actionToggleCB ActionButton CBrite}; % ... find the ActionButton method actionToggleCB
else, % simple action button, no toggle
end

hbut = uicontrol(hp, 'units', 'pixels', 'position', [X Y Q.Extent], ...
    'style', 'pushbutton', ...
    'BackgroundColor', get(hp,'BackgroundColor'), ...
    'Callback', CBleft, ...
    'ButtonDownFcn', CBrite, ...
    'ForegroundColor', get(hp,'ForegroundColor'), ...
    'String', Q.CurrentString, 'TooltipString', Q.Tooltip, ...
    Q.uicontrolProps);

figh = parentfigh(hp);
figh = figh.Number;
% add menuitem with accelerator key if needed
if ~isempty(Q.Accel),
    mh = findobj(figh, 'type', 'uimenu', 'label', Q.Accel.MenuLabel);
    if isempty(mh), %provide uimenu
        mh = uimenu(figh, 'label', Q.Accel.MenuLabel);
    end
    hmenu=uimenu(mh,'label', Q.Name, 'accelerator', Q.Accel.Key, ...
        'enable', get(hbut,'enable'), 'callback', {Q.Callback 'Left'});
end

% add Action info (including handles!) to userdata of GUI figure
Q.uiHandles.Button = hbut;
if ~isempty(Q.Accel), Q.uiHandles.Menu = hmenu; end;
Button = getGUIdata(figh, 'ActionButton', []);
Button = [Button Q];
setGUIdata(figh, 'ActionButton', Button);
% also set button's own userdata to Q
set(hbut,'userdata',Q);
if ~isempty(Q.Accel), set(hmenu,'userdata',Q); end;


%               Name: 'StartFrequency'
%             Prompt: 'Start:'
%             String: '12000 12001'
%               Unit: 'Hz'
%         Constraint: 'rreal/positive'
%            Tooltip: 'First frequency of series.'
%             MaxNum: 2
%     uicontrolProps: [1x1 struct]
% % % % % %         Units: 'pixels'
% % % % % %      FontSize: 8
% % % % % %      FontName: 'MS Sans Serif'
% % % % % %     FontUnits: 'points'




