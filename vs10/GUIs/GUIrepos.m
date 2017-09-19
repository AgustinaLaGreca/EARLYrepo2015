function GUIrepos(h);
% GUIrepos - reposition GUI window 
%    GUIrepos(h) restores the position of the GUI to the value it had when
%    calling GUIclose.
%
%    See also GUIclose, gcg.

if ~isGUI(h), error('handle does not belong to a GUI.'); end

% prepare saving of position @ close time
Name=getGUIdata(h,'GUIname');
placefig(h, Name);


