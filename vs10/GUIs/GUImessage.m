function GUImessage(figh, Mess, Mode, Qnames);
% GUImessage - display general message in GUI
%    GUImessage(figh, Mess, Mode) displays text message Mess in the GUI
%    with handle figh. The text will be displayed in the Main Messenger 
%    of the GUI. Mess may be a char string or cellstring. Mode is one
%    of 'shy', 'neutral', 'warning', 'error', 'append'. The Mode affects the 
%    color of the text displayed as prescribed by GUIsettings. The display
%    of 'shy' is the same as that of 'neutral', but for 'shy' the GUI is
%    not brought to the fore. This is useful when calling GUImessage from a
%    timer or other background process: when using the other modes, the 
%    commandline disappears on each broadcast. The default mode is 'neutral'. 
%    The 'append' mode appends the text to the text already displayed. 
%
%    If Mess is empty, nothing is done, and any text currently displayed 
%    is unaffected. Specify a single blank, ' ', to empty the display.
%
%    GUImessage(figh, Mess, Mode, Qnames) also sets the background color of 
%    the edit control(s) and/or toggle buttons belonging to queries om figh
%    whose names are in cellstring Qnames. Their color is set to the 
%    "error color" (see paramquery/coloredit). If Mess is empty, he edit 
%    control is left alone. Any Postfixes 'Edit' are stripped off the names.
%
%    If figh is not an existing figure handle, the message is simply
%    displayed on the command screen.
%
%    See also Messenger, paramquery/coloredit, GUIsettings.

if nargin<3, Mode = 'neutral'; end
if nargin<4, Qnames={}; end

if ~isSingleHandle(figh),
    if isequal('neutral', Mode), Mode = 'message'; end
    disp(['-----' Mode '-----']);
    if iscell(Mess), Mess = strvcat(Mess{:}); end
    disp(Mess);
    disp(' ');
    return
end

if isempty(Mess), return; end

% blink & redden edit controls.
if ischar(Qnames), Qnames = {Qnames}; end  % make it a cellstring
if ~iscellstr(Qnames),
    error('Input arg Qnames must be cell string with names of edits and toggles.');
end
Qnames = strrep(Qnames,'Edit', '');
Qnames = strrep(Qnames,'Unit', '');
Q = getGUIdata(figh,'Query', ParamQuery());
Q = Q(Qnames{:});
coloredit(Q,false); 
colortoggle(Q,false); 

persistent GSM; if isempty(GSM), GSM = GUIsettings('Messenger'); end

[Mode, ModeMess] = keywordMatch(Mode,{'neutral', 'warning', 'error', 'append' 'shy'}, 'Message mode ');
error(ModeMess);

M = getGUIdata(figh,'MainMessenger');
htext = M.uiHandles.Text; 
OldText = get(htext,'String'); OldText = OldText(:).'; % make sure OldText is a row array
Color = get(htext,'ForegroundColor');
if ~iscellstr(Mess), Mess = cellstr(Mess); end
Mess = Mess(:).'; % make sure Mess is a row array

DoFocus = 1;
switch Mode,
    case 'append',
        Mess = [OldText Mess];
    case {'neutral','shy'},
        Color = GSM.NeutralColor;
        DoFocus = ~isequal('shy', Mode) && ~isempty(MyFlag('neverfocus'));
    case 'warning',
        Color = GSM.WarnColor;
    case 'error',
        Color = GSM.ErrorColor;
end
set(htext,'string', Mess, 'foregroundcolor', Color);
drawnow; 
if DoFocus, figure(figh); end





