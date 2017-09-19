function M=CycleList(Name, Nmax, Callback, varargin);
% CycleList - constructor for CycleList objects.
%    L=CycleList(Name, UImenuLabel, Nmax, Callback, ...) creates, but does not 
%    render, a CycleList object, which can be added to a PulldownMenu.
%    Nmax is the maximum number of menu items in the list. When
%    adding an item (see CycleList/additem), it is placed on top of the
%    list, moving the existing items one down the list. The (Nmax+1)-th 
%    item, if any, is thereby removed. Callback is the callback function of
%    the menuitems to be added. Callback must be a single function handle
%    or function name. Cell arrays are not allowed! When a menu item of L
%    is clicked, the Callback function is called with four arguments:
%         feval(Src, Event, L.Callback, L, Item)
%    where Src is the handle of the calling uimenu item; Event is that void
%    thing that will mean something in the Matlab future; L is the actual, 
%    fully up-to-date version of the CycleList object L at the time of the 
%    callback; and Item is that element of L.Items that was selected by the 
%    mouse click. See StimGuiFilePulldown for a sample implementation.
%    The "..." are property/value pairs or struct containing uimenu 
%    properties used for rendering the uimenu items in the GUI.
%
%    CycleLists cannot be placed directly in a GUIpiece, but have to be
%    embedded in a PulldownMenu object, using the PulldownMenu/additem
%    method.
%
%    See also cycleList/additem PulldownMenu/additem.

%-----check input args-------

% specials cases: struct input or object input
if nargin==1 && isstruct(Name),
    M = Name;
    M = class(M, mfilename);
    return;
elseif nargin==1 && isa(Name, mfilename),
    M = Name;
    return;
elseif nargin<1, % void object with correct fields
    [Name, Nmax, Handle, parentHandle, ItemHandles, Items, Callback, uimenuProps] = deal([]);
    M = CollectInStruct(Name, Nmax, Handle, parentHandle, ItemHandles, Items, Callback, uimenuProps);
    M = class(M, mfilename);
    return;
end

%------regular call from here: AxesDisplay is fully specified ------------
if ~isvarname(Name),
    error('Name of cycleList object must be valid Matlab variable name (see ISVARNAME).');
end
error(numericTest(Nmax,'posint','Input argument Nmax is'));
if ~isfhandle(Callback) && ~isvarname(Callback),
    error('Callback of cycleList object must be function name or function handle. Cell arrays are not allowed.');
end

GSC = GUIsettings('CycleList');
% select those fields of GSP and varargin that match uimenu properties
DefaultProps = FullFieldnames(GSC,uimenuProperties, 'select');
ExplicitProps = FullFieldnames(struct(varargin{:}),uimenuProperties);
% merge default & explicit uipanel properties. The latter take precedence
uimenuProps = structJoin(DefaultProps, ExplicitProps);

% provide empty initializations for rendering-related fields & create object
[Handle, parentHandle, ItemHandles, Items] = deal([]);
M = CollectInStruct(Name, Nmax, Handle, parentHandle, ItemHandles, Items, Callback, uimenuProps);
M = class(M, mfilename);




    
    