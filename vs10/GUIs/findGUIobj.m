function M=findGUIobj(figh, Type, Name);
% findGUIobj - find named messenger object in GUI
%    findGUIobj(figh, 'Type', 'Name') returns the rendered GUI object
%    of specified Type named 'Name' if it exists in the GUI with handle figh. 
%    An error occurs if Foo does not exist or if object of the Type do not 
%    have a Name property. Type and Name are case INsensitive and may be 
%    abbreviated.
%
%    Use getGUIdata to find out what types of GUI objects exist in a GUI.
%
%    Examples:
%      findGUIobj(gcg, 'Query', 'SPL')
%      findGUIobj(gcg, 'Messenger', 'BulletinBoard')
%
%    See also GCG, messenger/report, GUIval.

gd = getGUIdata(figh);
[Type, Mess] = keywordMatch(Type, fieldnames(gd), 'GUIobject Type');
error(Mess);
M = gd.(Type);
MS = struct(M);
if ~isfield(MS, 'Name'), 
    error(['GUIobjects of type ''' Type ''' do not have a Name property.']);
end
[Name, Mess] = keywordMatch(Name, {MS.Name}, 'GUIobject Name');
error(Mess);
M = M(Name);




