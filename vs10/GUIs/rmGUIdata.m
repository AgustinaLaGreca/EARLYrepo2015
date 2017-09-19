function rmGUIdata(h, varargin);
% rmGUIdata - remove struct field of userdata property of a graphics object
%    rmUserdata(h, FN) removes the FN field of struct-valued 
%    'userdata' property. If the value of FN is a anction object, its
%    address is removed from the GUI's action list.
%
%    rmGUIdata(h) sets the entire userdata of h to []; If h has a 'tag' 
%    property set to 'has_userdata', it will also be set to [].
%
%    rmGUIdata([h1 h2 ..], ..) visits all handles one-by-one.
%
%    See also getGUIdata, GUIactionList, setGUIdata, findobj.

if length(h)>1, % recursive handling
    for ii=1:length(h),
        rmUserdata(h(ii), varargin{:});
    end
    return
end
%----------single h from here------------

error(handleTest(h,'any', 'h argument'));
ud = get(h,'userdata');
if ~isempty(ud) & ~isstruct(ud),
    error('Userdata contains non-struct data.');
end
if isstruct(ud) & numel(ud)>1,
    error('Userdata contains struct array with multiple elements.')
end

% special case: remove all GUIdata
if nargin<2, % delete all
    set(h,'userdata', []);
    if isequal('has_userdata', get(h,'tag'))
        set(h,'tag','');
    end
    return;
end

% specic field or subfield from here
FN = varargin{1};
% check for staggered fields
Mult_FN = Words2cell(FN,'.');
Staggered = numel(Mult_FN)>1;
if Staggered,
    nolast = cell2words(Mult_FN(1:end-1),'.'); % FN w/o last .field
    y = ud.(nolast); % first extracting ud.foo.fld.. gives overloaded subsref a chance to kick in. Direct use ...
    y = rmfield(y, Mult_FN{end}); % ... rmfield would bypass these overloaded subsref cases
    ud.(nolast) = y;
else, % remove single field
    if ~isfield(ud, FN), error(['Userdata field ''' FN ''' not present.']); end
    ud = rmfield(ud,FN);
end
set(h,'userdata',ud);
GUIactionList(h, FN, 'remove'); % this is faster than having to check whether FN's former value was indeed an action object





