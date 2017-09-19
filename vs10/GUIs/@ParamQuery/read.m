function [Value, Unit, Mess, hh]=read(Q, Flag);
% ParamQuery/read - read numerical value and unit from rendered paramquery.
%   [Value, Unit, Mess, h]=read(Q) reads the Value and Unit from the
%   edit and unit uicontrols of the paramaquery Q drawn in the GUI.
%   When the text typed in the edit field does not meet the constraints of
%   Q, Mess contains an nonempty error message, and Value and Unit are
%   empty. The last output arg, h, is a handle to the edit uicontrol,
%   which may be used to highlight any faulty user input. For edit-less
%   paramqueries, h is the handle of the unit button.
%
%   Special constraint values are
%       'rseed': a valid input to SetRandState (nan or positive int<Nmax)
%     'varname': a valid Matlab variable name (see ISVARNAME). 
%       'anwin': a pair of nonnegative numbers or the string 'burstdur'
%      'string': any char string
%   'togglenum': the number displayed on the toggle button
%
%   End users may suppy their own reading method Foo by choosing a
%   constraint value ~Foo when creating the paramquery. In this case, the 
%   user has to supply a paramquery method read_foo with the following
%   signature: [Value, EditStr, Mess]=local_varname(hh, Q). Value is the
%   return Value; EditStr is the (trimmed) string as read from the edit
%   field; Mess is empty unless something went wrong; then it contains a
%   message to be displayed; hh is the graphics handle to the edit
%   uicontrol; Q is the paramquery object as passed to paramquery/read.
%
%   Numerical values in the edit control may be specified as "k numbers",
%   e.g., 2k3 == 2300. On successful reading from an edit control, any
%   k-string is replaced by its full representation.
%
%   [Value, Unit, Mess]=read(Q, 'coloredit') sets the background color of
%   the edit uicontrol according to the result of the read action. The
%   colors used are determined by the following fields of the ParamQuery 
%   GUIsettings:
%        EditOkayColor: background color of edit uicontrol in case of success
%        EditErrorColor: background color of edit uicontrol in case of problems
%   For toggles, the color of the pushbutton is reset to its default
%   color.
%
%   See also GUIsettings, ParamQuery, ParamQuery/grab.

if nargin<2, Flag = '---'; end
[Flag, Mess] = keywordMatch(Flag,{'coloredit' '---'}, 'flag input argument');


if isempty(Q.uiHandles) || ~isfield(Q.uiHandles, 'Unit') ...
        || ~ishandle(Q.uiHandles.Unit),
    error(['Paramquery ''' Q.Name ''' is not currently rendered in a GUI.']);
end

Mess = ''; % optimistic default

% unit is simply the string prop of unit uicontrol
Unit = get(Q.uiHandles.Unit, 'string');
if isequal('coloredit', Flag),
    if istoggle(Q), 
        colortoggle(Q,1);
    end
end

% get value from edit string, if present
hh = getFieldOrDefault(Q.uiHandles, 'Edit', []);
SpecialSyntax = 0;
if isempty(hh) && ~isequal('togglenum', Q.Constraint), % no edit uicontrol: return the string of the unit toggle button
    Value = Unit;
    Unit = '';
    hh = getFieldOrDefault(Q.uiHandles, 'Unit', []);
    return; % value from toggle button: always okay
end

% reading from an edit uicontrol
SpecialSyntax = 1; % assume an exception; they are dealt with first
SkipNumelTest = 0; % assume that # elements still needs testing
if isequal('rseed', Q.Constraint),
    [Value, EditStr, Mess]=local_rseed(hh, Q);
elseif isequal('varname', Q.Constraint),
    [Value, EditStr, Mess]=local_varname(hh, Q);
elseif isequal('anwin', Q.Constraint),
    [Value, EditStr, Mess]=local_anwin(hh, Q);
    SkipNumelTest = 1;
elseif isequal('string', Q.Constraint),
    [Value, EditStr, Mess]=local_string(hh, Q);
elseif isequal('togglenum', Q.Constraint),
    hh = Q.uiHandles.Unit;
    [Value, Unit, EditStr, Mess]=local_togglenum(hh, Q);
elseif isequal('~', Q.Constraint(1)), % user-supplied interpreter
    UserMethod = fhandle(Q.Constraint(2:end));
    [Value, EditStr, Mess]=feval(UserMethod, hh, Q);
    SkipNumelTest = 1; % assume user defined function works properly
else, % regular case: read a number or multiple numbers
    SpecialSyntax = 0;
    [Value, EditStr] = kstr2num(get(hh,'string'));
    if any(isnan(Value)),
        Mess = ['Non-numerical value of ' Q.Name '.'];
    end
end
% further processing of values obtained from edit control

if isempty(Mess) && ~SpecialSyntax, % perform further tests
    Mess = numericTest(Value, Q.Constraint, [Q.Name ' is ']);
end

if isempty(Mess) && ~SkipNumelTest && (numel(Value)>Q.MaxNum),
    Mess = [Q.Name ' has too many (>' num2str(Q.MaxNum)  ') values.'];
end

if isempty(Mess),
    set(hh, 'String',  EditStr);
end

persistent GSQ; if isempty(GSQ), GSQ = GUIsettings('ParamQuery'); end

%        OkayColor: background color of edit uicontrol in case of success
%        ErrorColor: background color of edit uicontrol in case of problems

if isequal('coloredit', Flag),
    coloredit(Q,isempty(Mess)); % if not okay: blink & redden
end


%================================
function [Value, Str, Mess]=local_rseed(hh, Q);
Mess = '';
Str = trimspace(get(hh,'string'));
if isempty(Str),
    Mess = [Q.Name ' is not specified.']
    Value = [];
elseif isequal('nan', lower(Str)), % fresh rand value from clock
    Value = SetRandState;
elseif length(Str)>=3 && isequal('exp', lower(Str(1:3))), % derive sed from current experiment
    Tail = 15634*str2num(['0' Str(4:end)]);
    ExpName = upper(name(current(experiment)));
    if isempty(ExpName),
        Mess = 'No current experiment defined.';
        Value = [];
    else,
        Value = Tail + Str2Seed(strrep(ExpName(max(1,end-5):end),'0', 'Q'));
        Value = mod(Value, 61^5);
    end
else, % must be integer value 0<I<Nmax
    [Value,Str] = kstr2num(Str);
    if isnan(Value), 
        Mess=['Non-numerical value of ' Q.Name '.']; 
    else,
        Mess = numericTest(Value, 'nonnegint', [Q.Name ' is ']);
    end
    if isempty(Mess),
        [dum , Nmax] = SetRandState([]);
        if Value>Nmax, Mess = [Q.Name ' exceeds maxumim value of ' num2str(Nmax) '.']; end
    end
end

function [Value, Str, Mess]=local_varname(hh, Q);
% read MatLab identifier
Mess = ''; Value = '';
Str = trimspace(get(hh,'string'));
if isempty(Str),
    Mess = [Q.Name ' is not specified.'];
    Value = [];
else, % must be valid Matlab identifier
    if ~isvarname(Str), 
        Mess=['Invalid identifier specified for ' Q.Name '.']; 
    elseif ismember(lower(Str),{'default', 'factory', 'remove'}),
        Mess=['''default'', ''factory'', and ''remove'' are invalid values.'];
    else, % okay
        Value = Str;
    end
end

function [Value, Str, Mess]=local_string(hh, Q);
% read string - always valid except when empty
Mess = ''; Value = '';
Str = trimspace(get(hh,'string'));
if isempty(Str),
    Mess = [Q.Name ' is not specified.'];
    Value = [];
elseif ismember(lower(Str),{'default', 'factory', 'remove'}),
        Mess=['''default'', ''factory'', and ''remove'' are invalid values.'];
else, % okay
    Value = Str;
end

function [Value, Unit, Str, Mess]=local_togglenum(hh, Q);
% take string from unit button, remove unit and convert to number
Mess = ''; Value = '';
Str = get(hh,'string');
W = Words2cell(Str,' ');
Mess = '';
Value = str2num(W{1});
Unit = W{2};

function [Value, Str, Mess]=local_anwin(hh, Q);
% read analysis window: pair of numbers or "burstdur"
Mess = ''; Value = '';
Str = trimspace(get(hh,'string'));
if isempty(Str),
    Mess = [Q.Name ' is not specified.'];
    Value = [];
elseif ~isempty(strmatch(lower(Str), 'burstdur')), 
    Value = 'burstdur';
else, %pair of values
    Value = str2num(Str);
    if numel(Value)~=2,
        Mess = [Q.Name ' must be pair of numbers or the string ''burstdur''.'];
        Value = [];
    elseif any(Value<0),
        Mess = [Q.Name ' elements may not be negative'];
        Value = [];
    elseif diff(Value)<=0,
        Mess = [Q.Name ' must be finite interval t0 t1 with t0<t1.'];
        Value = [];
    end
end


