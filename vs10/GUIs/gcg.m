function h=gcg;
% GCG - handle to current GUI
%    GCG returns the handle to the most recently used GUI, i.e., that GUI
%    which has been used most recently (the one on top of the stack).
%    A GUI is a figure created by newGUI. If no GUI is open, gcg returns
%    [].
%
%    See also newGUI, gcf, issinglehandle.


h = [];
% get all figure handles, also hidden ones
shh = get(0,'showhiddenhand'); % store to restore
set(0,'showhiddenhand', 'on');
hf = findobj(0,'type', 'figure');
hf = hf.Number;
set(0,'showhiddenhand', shh);

for ii=1:numel(hf),
   nam = getGUIdata(hf(ii), 'GUIname',nan);
   if ischar(nam), % gotcha
       h = hf(ii);
       break;
   end
end







