function Q=draw(hp, Q, XY);
% AxesDisplay/draw - draw AxesDisplay object.
%   draw(hp, B, X,Y) draws AxesDisplay object B in uipanel with handle hp at
%   position (X,Y) in pixels. This function is typically recursively called 
%   by GUIpiece/draw.
%
%   See AxesDisplay, GUIpiece/draw, GUIpiece.

cpos=get(hp,'position');
X=XY(1); Y=cpos(4)-XY(2)-Q.Extent(2);

hax = axes('parent', hp, 'units', 'pixels', 'outerposition', [X Y Q.Extent], ...
    Q.axesProps);


% add Action info (including handles!) to userdata of GUI figure
Q.uiHandles.Axes = hax;
figh = parentfigh(hp);
AD = getGUIdata(figh, 'AxesDisplay', []);
AD = [AD Q];
setGUIdata(figh, 'AxesDisplay', AD);
% also set button's own userdata to Q
set(hax,'userdata',Q);

if isequal('@', Q.Name(1)), % main AxesDisplay
    setGUIdata(figh, 'MainAxesDisplay', Q);
end



