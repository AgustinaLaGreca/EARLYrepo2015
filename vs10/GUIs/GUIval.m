function S=GUIval(figh,varargin);
% GUIval - evaluate user input from GUI
%    S=GUIval(h) reads all paramqueries displayed on the GUI with handle h
%    and returns the result in a struct S. For each paramquery named Foo,
%    two fields are created in S: 
%        Foo: numerical value of parameter Foo
%        FooUnit: string containing the unit of Foo (e.g. 'Hz')
%    Finally, the handles of all edit controls are returned in struct
%    S.handle. If Foo does not have an edit control, the handle to its
%    unit control is returned instead.
%
%    If a faulty value is encountered (e.g. empty edit field, 
%    non-numerical value), GUIval skips the remaining queries, and returns 
%    the value [] for S. Also an appropriate error message is displayed 
%    using GUImessage and the offending edit fields are highlighted. The 
%    highlighting is undone upon the next GUIval call.
%
%    See also GUIgrab, paramquery/read, GUImessage, paramquery/highlight.


figure(figh);
S.GUIname = getGUIdata(figh,'GUIname');
S.GUIhandle = figh;
S.CreatedBy = getGUIdata(figh,'CreatedBy');

Q = getGUIdata(figh,'Query'); % array of queries.
% read queries one by one. Stop in case of problems.
Handle.GUIfig = figh;
Mess = ''; 
for ii=1:numel(Q),
    q = Q(ii);
    if ~isempty(separator(q));
        SeP = separator(q);
        S.(['i_________' SeP '_______']) = ['_____' SeP '_____'];
    end
    [Value, Unit, Mess h]=read(q,'coloredit');
    if ~isempty(Mess),
        if nargin<2
            break;
        end
    end
    S.(q.Name) = Value;
    S.([q.Name 'Unit']) = Unit;
        Handle.(q.Name) = h;
end
S.handle = Handle;
if isempty(Mess), 
    S.GUIgrab = GUIgrab(figh); % store unprocessed user input for reference
    S.z_______________ = '_________________'; % separator
else, % report error & return empty S
    GUImessage(figh,Mess,'error');
    S = []; 
end



