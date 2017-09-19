function Q=draw(hp, Q, XY);
% Messenger/draw - draw Messenger object.
%   draw(hp, B, X,Y) draws Messenger object B in uipanel with handle hp at
%   position (X,Y) in pixels. This function is typically recursively called 
%   by GUIpiece/draw.
%
%   See Messenger, GUIpiece/draw, GUIpiece.

cpos=get(hp,'position');
X=XY(1); Y=cpos(4)-XY(2)-Q.Extent(2);

htext = uicontrol(hp, 'units', 'pixels', 'position', [X Y Q.Extent], ...
    'style', 'text', ...
    'BackgroundColor', get(hp,'BackgroundColor'), ...
    'String', Q.TestLine, ...
    'HorizontalAlignment', 'left', ...
    Q.uicontrolProps);


% add Action info (including handles!) to userdata of GUI figure
Q.uiHandles.Text = htext;
figh = parentfigh(hp);
figh = figh.Number;
Messenger = getGUIdata(figh, 'Messenger', []);
Messenger = [Messenger Q];
setGUIdata(figh, 'Messenger', Messenger);
% also set button's own userdata to Q
set(htext,'userdata',Q);

if isequal('@', Q.Name(1)), % main messager
    setGUIdata(figh, 'MainMessenger', Q);
end


