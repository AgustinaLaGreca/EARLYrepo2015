function P = additem(P, Label, Callback, varargin);
% PulldownMenu/additem - add item to PulldownMenu
%   P=additem(P, Label, Callback, ...) adds a menuitem to PulldownMenu P.
%   Label and Callback refer to the Label and Callback properties of the
%   uimenu realizing the item. "..." stands for any number of
%   property/value pairs (or a property/value struct) to be imposed on the
%   uimenu items. Callback may be a string, a function handle, or a cell
%   array containing a function handle and addition arguments (See Matlab
%   documentation on callbacks).
%
%   P=additem(P, p), where p is a PuldownMenu object, adds pulldown menu p
%   to P. When rendered, p will show when its item in P is clicked.
%
%   P=additem(P, CL), where CL is a Cyclelist object, adds a cycle list to
%   the menu.
%
%   Note: Items must be added to P *before* rendering P. A single draw(h,P)
%   will then recursively draw all items contained in P.
%
%   See also PulldownMenu, CycleList, PulldownMenu/draw.


if nargout<1, error('Too few output args'); end
if isSingleHandle(P.Handle),
    error('Cannot add items to pulldownMenu object that is already rendered');
end
if isvoid(P),
    error('Cannot add items to void PulldownMenu.');
end

% add the item to the list in P. Distinguish the 3 cases (see help text)
switch lower(class(Label)),
    case 'char', % "elementary" item
        newItem = CollectInStruct(Label, Callback);
        if ~isempty(varargin),
            XProps = struct(varargin{:});
            XProps = FullFieldnames(XProps, uimenuProperties, 'select');
            newItem = structJoin(newItem, XProps);
        end
    case {'pulldownmenu' ,'cyclelist'}, % recursive pulldown menu
        newItem = Label;
        if nargin>2, error('Too many input args.'); end
    otherwise,
        error('Second input arg must be string (Label), or PulldownMenu or CycleList object.');
end
P.Items{end+1} = newItem; % strangely, P.Items = [P.Items newItem] fails








