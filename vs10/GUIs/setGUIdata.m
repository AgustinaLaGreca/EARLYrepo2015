function P=setGUIdata(h, FN, Value, flag);
% setGUIdata - set struct field of userdata property of a graphics object
%    setGUIdata(h, FN, Value) "uploads" a Value to a GUI under the name of
%    FN. That is, it sets the FN field of struct-valued 
%    'userdata' property of h to Value. If h has a 'tag' property, 
%    it will be set to 'has_userdata'. setGUIdata returns a GUIaccess
%    struct that gives access to the FN variable in the way described in
%    GUIaccess.
%
%    The names FN may also contain subfields, like 'Foo.fld', etc, in 
%    which case the 'fld' subfield of Foo is set. 
%
%    A=setGUIdata(h, FN, Value) returns the 'address' A of FN within the GUI
%    as described in GUIaccess.
%
%    If Value is an action object, setGUIdata(figh, FN, Value) has the side
%    effect that its address is stored in a list of action objects uploaded 
%    to the GUI. See GUIactionList for details.
%    A=setGUIdata(h, FN, Value, 'dontlist') prevents FN from being listed
%    this way.
%
%    See also getGUIdata, rmGUIdata, GUIaccess, GUIactionList.

if nargin<4, flag=''; end

error(handleTest(h,'any', 'h argument'));

% if Value is an Action object, initialze action list here, before FN itself so that in getGUIdata(figh) ..
if ~isequal('dontlist', flag) && isa(Value, 'action'), % .. the ActionList__ item appears above its members.
       GUIactionList(h, 'init');
end
ud = get(h,'userdata');
if ~isempty(ud) && ~isstruct(ud),
    error('Userdata contains non-struct data.');
end
if isstruct(ud) && numel(ud)>1,
    error('Userdata contains struct array with multiple elements.')
end

% check for staggered fields
Mult_FN = Words2cell(FN,'.');
Staggered = numel(Mult_FN)>1;

if isempty(ud),
    clear ud
    ud.(FN) = Value;
elseif Staggered, % struct.struct..
    eval(['ud.' FN ' = Value;']); % the old way of setting struct fields
else, % struct
    ud(1).(FN) = Value;
end
set(h,'userdata', ud);
P = GUIaccess(h,FN);

try,
    set(h,'tag', 'has_userdata');
end % try

if ~isequal('dontlist', flag) && isa(Value, 'action'), % add to the action list
       GUIactionList(h, FN);
end
