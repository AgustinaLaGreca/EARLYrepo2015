function GUIclose(h);
% GUIclose - close GUI window after saving its position
%    GUIclose(h) first checks whether the GUI with handle h is closable 
%    (see GUIclosable). If not, nothing is done. If the GUI is closable,
%    the GUI is closed after saving its position.
%    Next time, the position can be restored by GUIrepos.
%
%    See also GUIclosable, GUIrepos, gcg.

if ~GUIclosable(h), return; end;
qq=getGUIdata(h,'SavePosAndDelete');
feval(qq{:});



