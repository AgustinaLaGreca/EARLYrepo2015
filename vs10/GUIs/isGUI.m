function ig = isGUI(h);
% isGUI - true for handle of a GUI
%    isGUI returns true (1) if h is the handle of a GUI.
%    In order to be a GUI handle, h has to be a figure handle, and the
%    GUIdata in h should contain a GUIname field.
%
%    See also newGUI, gcg.

ig = ishandle(h) && isequal(get(h,'type'), 'figure') ...
    &&     ~isempty(getGUIdata(h,'GUIname',[]));




