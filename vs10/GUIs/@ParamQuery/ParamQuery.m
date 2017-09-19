function P=ParamQuery(Name, Prompt, String, Unit, Constraint, Tooltip, MaxNum, varargin);
% ParamQuery - constructor for ParamQuery objects.
%    P=ParamQuery(Name, Prompt, String, Unit, Constraint, Tooltip, MaxNum, ...);
%    "..." are property/value pairs or struct containing uicontrol
%    properties used for rendering the query in a GUI.
%    Inputs are
%        Name: varname string. This name is used in error messages when
%              evaluating the GUI (see GUIval).
%      Prompt: char string displayed to the left of the edit control.
%      String: character string determining the size (extent) of the edit
%              control. If String is empty, no edit uictrol will be
%              displayed.
%        Unit: physical unit of the parameter. Will be displayed to the
%              right of the edit control. If String is cell string instead
%              of a char string, the units will be rendered as a toggle
%              button, leaving the choice of units to the user.
%    Constraint: value constraints that are tested when reading the value
%             of rendered paramquery objects. (see paramquery/read). The
%             string is either as described in NumericTest, or a one of a
%             fixed set of special values  rseed|varname, or a string of
%             the form '~foo'. In the latter case, reading will be
%             delagated to a user-supplied paramquery method
%             paramquery/read_foo.
%     Tooltip: tooltip string shown at the prompt of the rendered query.
%             Cell strings are converted to char(10)-delimited strings.
%      MaxNum: maximum number of values that can be specified by the user.
%              Defaults to one.
%
%    Examples:
%    Nreps = paramquery('Nreps', '# Reps:', '1500', '', ...
%            'rreal/posint', 'Number of repetitions of each condition.',1);
%    RSeed = paramquery('RSeed', 'Rand Seed:', '1234567', '', ...
%      'rseed', {'Random seed used for presentation order.' 'Specify NaN to refresh seed upon each usage.'},1);
%
%    See also GUIval, paramquery/read.


%-----check input args-------

% specials cases: struct input or object input
if nargin==1 && isstruct(Name),
    P = Name;
    P = class(P, mfilename);
    return;
elseif nargin==1 && isa(Name, mfilename),
    P = Name;
    return;
elseif nargin<1, % void object with correct fields
    [Name, Extent, Prompt, String, Unit, Constraint, Tooltip, MaxNum, uiHandles, uicontrolProps, Parent, Separator] = deal([]);
    P = CollectInStruct(Name, Extent, Prompt, String, Unit, Constraint, Tooltip, MaxNum, uiHandles, uicontrolProps, Parent, Separator);
    P = class(P, mfilename);
    return;
end

%------regular call from here: Param is fully specified except maybe MaxNum------------
if nargin<7, MaxNum=1; end
if ~isvarname(Name),
    error('Name of paramquery object must be valid Matlab variable name (see ISVARNAME).');
end
if ~ischar(Prompt),
    error('Prompt property of paramquery object must be character string.');
end
if ~ischar(String),
    error('String property of paramquery object must be character string.');
end
if iscellstr(Tooltip),
    Tooltip = cellstr2char10str(Tooltip);
end
if ~ischar(Tooltip),
    error('Tooltip property of paramquery object must be cellstring or character string.');
end
if ~ischar(Unit) && ~iscellstr(Unit),
    error('Unit property of paramquery object must be character string or cellstring.');
end
if ~ischar(Constraint),
    error('Constraint property of paramquery object must be character string.');
end
if ~isnumeric(MaxNum) || ~isscalar(MaxNum) || ~isequal(round(MaxNum), MaxNum) || (MaxNum<=0),
    error('MaxNum property of paramquery object must be a positive integer.');
end

persistent GSP; if isempty(GSP), GSP = GUIsettings('Paramquery'); end
DefaultParamqueryProps = FullFieldnames(GSP,uicontrolProperties,'select');
ExplicitParamqueryProps = FullFieldnames(struct(varargin{:}),uicontrolProperties);
uicontrolProps = combineStruct(DefaultParamqueryProps, ExplicitParamqueryProps);
Extent = [];
uiHandles = [];
Parent = [];
Separator = '';
P = CollectInStruct(Name, Extent, Prompt, String, Unit, Constraint, Tooltip, ...
    MaxNum, uiHandles, uicontrolProps, Parent, Separator);
P = class(P, mfilename);

P.Extent = [width(P) height(P)];



    
    