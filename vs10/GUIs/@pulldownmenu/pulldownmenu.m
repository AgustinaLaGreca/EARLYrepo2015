function M=PulldownMenu(Name, Label, varargin);
% PulldownMenu - constructor for PulldownMenu objects.
%    P=PulldownMenu(Name, UImenuLabel, ...) creates, but does not 
%    render, a PulldownMenu object, which can be added to a GUI.
%    UImenuLabel is the label of the uimenu in the GUI. 
%    The "..." stand for any number of uimenu property/value pairs or a
%    struct to be used when creating the items under P.
%
%    Use the method PulldownMenu/additem to add uimenu items to P.
%
%    See also PulldownMenu/additem, "methods PulldownMenu".


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
    [Name, Label, Handle, ItemHandles, CycleItems, Items, uimenuProps] = deal([]);
    M = CollectInStruct(Name, Label, Handle, ItemHandles, Items, CycleItems, uimenuProps);
    M = class(M, mfilename);
    return;
end

%------regular call from here: AxesDisplay is fully specified ------------
if ~isvarname(Name),
    error('Name of Pulldown object must be valid Matlab variable name (see ISVARNAME).');
end
if ~ischar(Label),
    error('Input argument uimenuLabel must be char string.');
end

GSP = GUIsettings('PulldownMenu');
% select those fields of GSP and varargin that match uimenu properties
DefaultProps = FullFieldnames(GSP,uimenuProperties, 'select');
ExplicitProps = FullFieldnames(struct(varargin{:}),uimenuProperties);
% merge default & explicit uipanel properties. The latter take precedence
uimenuProps = structJoin(DefaultProps, ExplicitProps);

% provide empty initializations for rendering-related fields & create object
[Handle, ItemHandles, CycleItems] = deal([]);
Items = {};
M = CollectInStruct(Name, Label, Handle, ItemHandles, Items, CycleItems, uimenuProps);
M = class(M, mfilename);




    
    