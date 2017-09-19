function ic=GUIclosable(h,ic);
% GUIclosable - get/set closable state of GUI
%    GUIclosable(h) returns the current closable state of GUI with handle h.
%    At creation time (see newGUI), the closable state is true (1).
%
%    GUIclosable(h,1) makes the GUI closable.
%
%    GUIclosable(h,0) makes the GUI unclosable.
%
%    See also GUIclose, GUIrepos, gcg.

if ~isGUI(h), error('handle does not belong to a GUI.'); end

if nargin<2, % get
    ic = getGUIdata(gcg,'OkayToClose');
else,
    setGUIdata(h,'OkayToClose',ic);
end



