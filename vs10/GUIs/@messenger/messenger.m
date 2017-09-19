function M=Messenger(Name, TestLine, Nlines, varargin);
% Messenger - constructor for Messenger objects.
%    M=Messenger(Name, TestLine, Nlines, ...)
%    TestLine is the initial string or cellstring, which determines the
%    extent of M. Nlines is the number of lines to be displayed.
%    "..." are property/value pairs or struct containing uicontrol 
%    properties used for rendering the text-style uicontrol in the GUI.
%
%    A special case is the "main messenger" to which all general messages
%    of a GUI are directed. To indicate that M is destined to be the main
%    messenger of the GUI, preprend its name with a '@' as in
%    M=Messenger(Name, TestLine, Nlines, FormatString, ...)


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
    [Name, Extent, TestLine, Nlines, uiHandles, uicontrolProps] = deal([]);
    M = CollectInStruct(Name, Extent, TestLine, Nlines, uiHandles, uicontrolProps);
    M = class(M, mfilename);
    return;
end

%------regular call from here: Messenger is fully specified ------------
IsMainMess = isequal('@', Name(1));
if ~isvarname(Name) && (IsMainMess && ~isvarname(Name(2:end))),
    error('Name of Messenger object must be valid Matlab variable name (see ISVARNAME).');
end
error(numericTest(Nlines, 'posint','Nlines input argument is '));
if ~ischar(TestLine) && (Nlines>1),
    error('If Nlines>1, TestLine property of Messenger object must be character string.');
elseif ~ischar(TestLine) && ~iscellstr(TestLine),
    error('TestLine property of Messenger object must be cellstring or character string.');
end

GSM = GUIsettings('Messenger');
% select those fields of GSP and varargin that match uicontrol properties
DefaultUicontrolProps = FullFieldnames(GSM,uicontrolProperties, 'select');
ExplicitUicontrolProps = FullFieldnames(struct(varargin{:}),uicontrolProperties);
% merge default & explicit uipanel properties. The latter take precedence
uicontrolProps = structJoin(DefaultUicontrolProps, ExplicitUicontrolProps);

% SizeString determines size of edit control
%uicontrolProps
Extent = StringExtent(repmat(TestLine,Nlines,1), uicontrolProps);
TwoSpace = StringExtent('  ', uicontrolProps);
Extent(1) = Extent(1) + TwoSpace(1);
uiHandles = [];
M = CollectInStruct(Name, Extent, TestLine, Nlines, uiHandles, uicontrolProps);
M = class(M, mfilename);




    
    