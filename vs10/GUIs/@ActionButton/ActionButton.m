function B=ActionButton(Name, String, SizeString, Tooltip, Callback, varargin);
% ActionButton - constructor for ActionButton objects.
%    B=ActionButton(Name, String, SizeString, Tooltip, Callback, ...)
%    Name: name of B. Must be valid Matlab identifier (see ISVARNAME)
%    String: char string shown on button or cell array of strings. In
%        the latter case, clicking the button results in rotating the 
%        string of the button to the next element of the array. On
%        rendering the button, the first element is used. See
%        ActionToggleCB for details on toggling ActionButtons.
%    SizeString: virtual char string determining button size.
%       Tooltip: tooltip string of button.
%      Callback: function called when clicking the button. Callback may
%                either be a function handle or a char string. No cell 
%                arrays are allowed! In case of a function handle @Foo, Foo 
%                is called with three arguments: first two are handle and 
%                Event as provided by Matlab. Third arg is the string 'Left'
%                or 'Right' depending on whether the button was left- or 
%                right-clicked. Note that when the button is disabled, 
%                left-clicking results in calling the Foo with 3rd arg 
%                'Right'. Use the 'SelectionType' property of the figure to 
%                figure out what kind of click was actually performed. In 
%                case of a callback string 'Foo', the function Foo is 
%                called with no arguments. In case of a string 'Foo A B', 
%                Foo is called with char string arguments 'A' and 'B', etc. 
%                If you have to pass additional information to the callback
%                function, use B=setXdata(B,D) to associate data D with the
%                actionbutton. Within the callback function, these data are
%                retrieved by:
%                    B = get(gcbo, 'userdata'); % retrieve action button
%                    D = getXdata(B);
%    "...": property/value pairs or struct containing uicontrol
%       properties used for rendering the pushbutton in a GUI.
%
%    See also ActionButton/draw, ActionButton/setXdata, ActionButton/ActionToggleCB,
%             ActionButton/accelerator, ActionButton/highlight,
%             "methods ActionButton"

%-----check input args-------

% specials cases: struct input or object input
if nargin==1 && isstruct(Name),
    B = Name;
    B = class(B, mfilename);
    return;
elseif nargin==1 && isa(Name, mfilename),
    B = Name;
    return;
elseif nargin<1, % void object with correct fields
    [Name, Extent, String, SizeString, CurrentString, PreviousString, Tooltip, Callback, Accel, Xdata, uiHandles, uicontrolProps] = deal([]);
    B = CollectInStruct(Name, Extent, String, SizeString, CurrentString, PreviousString, Tooltip, Callback, Accel, Xdata, uiHandles, Extent, uicontrolProps);
    B = class(B, mfilename);
    return;
end

%------regular call from here: ActionButton is fully specified ------------
if ~isvarname(Name),
    error('Name of ActionButton object must be valid Matlab variable name (see ISVARNAME).');
end
if ~ischar(String) && ~iscellstr(String),
    error('String property of ActionButton object must be character string or cell array of strings.');
end
if ~ischar(Tooltip),
    error('Tooltip property of ActionButton object must be character string.');
end
if ~ischar(Callback) && ~isfhandle(Callback),
    error('Callback of actionButton must be function handle or char string. Cell arrays are not allowed.');
end

GSB = GUIsettings('ActionButton');
% select those fields of GSB and varargin that match uicontrol properties
DefaultUicontrolProps = FullFieldnames(GSB,uicontrolProperties, 'select');
ExplicitUicontrolProps = FullFieldnames(struct(varargin{:}),uicontrolProperties);
% merge default & explicit properties. The latter take precedence
uicontrolProps = structJoin(DefaultUicontrolProps, ExplicitUicontrolProps);

% SizeString determines size of button
Extent = StringExtent([SizeString '  '], uicontrolProps)*GSB.ExtentMatrix;
uiHandles = []; Accel = []; Xdata = [];
% initial value of string property
if iscellstr(String), CurrentString = String{1};
else, CurrentString = String;
end
PreviousString = '';
B = CollectInStruct(Name, Extent, String, SizeString, CurrentString, PreviousString, Tooltip, Callback, ...
    Accel, Xdata, uiHandles, ...
    Extent, uicontrolProps);
B = class(B, mfilename);




    
    