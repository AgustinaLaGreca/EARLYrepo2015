function M=AxesDisplay(Name, Extent, varargin);
% AxesDisplay - constructor for AxesDisplay objects.
%    M=AxesDisplay(Name, Extent, ...) creates, but does not draw, an
%    AxisDisplay object, which can be added to a GUIpiece.
%    Extent is the size (Width,Height) in pixels.
%    "..." are property/value pairs or struct containing axes
%    properties used for rendering the axes in the GUI.

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
    [Name, Extent, uiHandles, axesProps] = deal([]);
    M = CollectInStruct(Name, Extent, uiHandles, axesProps);
    M = class(M, mfilename);
    return;
end

%------regular call from here: AxesDisplay is fully specified ------------
IsMainAx = isequal('@', Name(1));
if ~isvarname(Name) && (IsMainAx && ~isvarname(Name(2:end))),
    error('Name of AxesDisplay object must be valid Matlab variable name (see ISVARNAME).');
end
error(numericTest(Extent, 'rreal/positive','Extent input argument is '));
if ~isequal(2,numel(Extent)),
    error('Extent input arg must be 2-element array.');
end

GSA = GUIsettings('AxesDisplay');
% select those fields of GSP and varargin that match uicontrol properties
DefaultProps = FullFieldnames(GSA,uicontrolProperties, 'select');
ExplicitProps = FullFieldnames(struct(varargin{:}),axesProperties);
% merge default & explicit uipanel properties. The latter take precedence
axesProps = structJoin(DefaultProps, ExplicitProps);

% SizeString determines size of edit control
%axesProps
uiHandles = [];
M = CollectInStruct(Name, Extent, uiHandles, axesProps);
M = class(M, mfilename);




    
    