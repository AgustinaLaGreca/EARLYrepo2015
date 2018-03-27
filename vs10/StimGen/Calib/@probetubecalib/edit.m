function [P, okay] = edit(P, kw, varargin);
% probetubecalib/edit - line editor for Probetubecalib objects.
%    PTC = edit(PTC) plots T and prompts the user for simple commands.
%    Commands are either a  single keyword or have the form  
%    "f1 f2 Action param". Examples:
%
%         4000 8000 Smooth 0.1
%    Smooth the 4000-8000 Hz region using a 0.1-Oct-wide window.
%
%         200 7000 Clip -40 120
%    Clip the 200-7000 Hz region between -40 dB and 120 dB.
%
%    After entering a valid command, a preview of the edited version is
%    shown and the user is asked to either accept or cancel the edit.
%    Note that the actions only apply to the "Tube" part of the
%    probetubecalib object.
%
%    After leaving the editor by pressing OK, the edited PTC is returned.
%    After leaving the editor by pressing CANCEL, the unaltered PTC is 
%    returned.
%
%    [PTC, Okay] = edit(PTC) also returns a logical Okay indicating whether
%    the OK button was pushed or not (i.e. whether PTC was changed or not).
%
%    Probetubecalib/Edit is also the callback function for the editors GUI.
%
%    See also Probetubecalib/smooth, Probetubecalib/clip.

if nargin<2, kw = 'launch'; end

GUIname = 'probetubecalibEDIT';
switch lower(kw),
    case 'launch', % edit(PTC, 'launch')
        figh = local_launch(P);
        GUIreaper(figh, 'wait'); % wait for permission to reap ...
        P = getGUIdata(figh, 'currentPTC'); % ... grab the data ...
        okay = getGUIdata(figh, 'wasEdited', false); %... and their status before ...
        local_delete(figh); % ... deleting the figure
    case 'keypress',
    case {'reset'}, % reset button
        LR = varargin{1};
        if ~isequal('Left', LR), return; end; % ignore right-clicks
        figh = gcbf;
        figh = figh.Number;
        Pi = getGUIdata(figh, 'initPTC');
        setGUIdata(figh, 'currentPTC', Pi);
        local_plot(figh);
    case {'submit'}, % submit button
        local_cmd(gcbf); 
    case {'ok'}, % OK button
        LR = varargin{1};
        if ~isequal('Left', LR), return; end; % ignore right-clicks
        figh = gcbf;
        figh = figh.Number;
        setGUIdata(figh, 'wasEdited', true); % notify that editing has been done
        GUIreaper(figh, 'reap'); % signal closing the figure; see 'launch' for data handling
    case {'cancel'}, % Cancel button
        LR = varargin{1};
        if ~isequal('Left', LR), return; end; % ignore right-clicks
        % replace current data by original data
        figh = gcbf;
        figh = figh.Number;
        Pi = getGUIdata(figh, 'initPTC');
        setGUIdata(figh, 'currentPTC', Pi);
        GUIreaper(figh, 'reap'); % signal closing the figure; see 'launch' for data handling
    case {'help'}, % editdum, 'Help', ...)
        local_help(varargin{:});
    case 'close', % same as cancel
        if nargin<3, figh=gcbf; else, figh=varargin{1}; end 
        figh = figh.Number;
        % close
        % replace current data by original data
        Pi = getGUIdata(figh, 'initPTC');
        setGUIdata(figh, 'currentPTC', Pi);
        GUIreaper(figh, 'reap'); % signal closing the figure; see 'launch' for data handling
    otherwise,
        error(['Unknown keyword ''' kw '''.']);
end

%=================================================================
%=================================================================
function figh=local_launch(PTC);
% launch ProbeTubeCalibration GUI
dum = probetubecalib(); % dummy Probetubecalib for recursive GUI calls
GUIname = 'probetubecalibEDIT';
CLR = 0.75*[1 1 1]+0.1*[0.4 -0.3 -0.6];
P_cmd = local_cmdPanel(CLR); % command
P_ax = local_actionPanel(CLR); % OK, CANCEL & messages
Hlp = pulldownmenu('Help', '&Help');
Hlp = additem(Hlp,'How to edit a probe tube calib data', @(Src,Ev)edit(dum, 'help', 'HowTo'));
%======GUI itself===========
PE=GUIpiece('ProbeTubeCalibEdit',[],[0 0],[10 4]);
PE = add(PE,P_cmd);
PE = add(PE, P_ax, below(P_cmd));
PE = marginalize(PE,[40 20]);
PE = add(PE, Hlp);
% open figure and draw GUI in it
figh = newGUI(GUIname, ['editing ' description(PTC)], {fhandle(mfilename), dum 'launch'}, 'color', CLR);
draw(figh, PE); 
% empty all edit fields & message field
GUImessage(figh, ' ','neutral');
% closereq & keypress fun
set(figh,'keypressf',{@(Src,Ev)edit(dum, 'keypress')});
set(figh,'closereq',{@(Src,Ev)edit(dum, 'close')});
% store PTC in GUIdata 
setGUIdata(figh,'initPTC', PTC);
setGUIdata(figh,'currentPTC', PTC);
% lauch plot and store handle
local_plot(figh);
% clear edits
Q = getGUIdata(figh, 'Query');
clearedit(Q);
%set(figh,'closereq')

%======
function S = local_cmdPanel(CLR); 
% panel for specifying stimulus params
dum = probetubecalib(); % dummy Probetubecalib for recursive GUI calls
S = GUIpanel('Stim', 'stimulus', 'backgroundcolor', CLR);
Cmd = ParamQuery('Cmd',  'Command:', '40000 45000 smooth 0.0025', '', 'string',  'Enter command such as "10e3 30e3 smooth 0.05". See help.', 1e3);
Submit = ActionButton('Submit', 'Submit', 'xxxxxxx', 'Apply command to probe tube trf.', @(Src,Ev,LR)edit(dum, 'Submit', LR));
S = add(S,Cmd);
S = add(S, Submit, nextto(Cmd));
S = marginalize(S,[3 5]);

%======
function Act = local_actionPanel(CLR);
% OK/Cancel & messages
dum = probetubecalib(); % dummy Probetubecalib for recursive GUI calls
Act = GUIpanel('Act', '', 'backgroundcolor', CLR);
MessBox = messenger('@MessBox', 'The problem is what you think it is or not ...',5, ... % the '@' in the name indicates that ...
    'fontsize', 12, 'fontweight', 'bold'); % MessBox will be the Main Messenger of he GUI
OK = ActionButton('OK', 'OK', 'XXXXXXXX', 'Accept edited version and exit.', @(Src,Ev,LR)edit(dum, 'OK', LR));
Cancel = ActionButton('Cancel', 'CANCEL', 'XXXXXXXX', 'Ignore all edits and quit.', @(Src,Ev,LR)edit(dum, 'Cancel', LR));
Reset = ActionButton('Reset', 'RESET', 'XXXXXXXX', 'Re-start with original version.', @(Src,Ev,LR)edit(dum, 'Reset', LR));
Act = add(Act,MessBox);
%==fast keys
OK = accelerator(OK,'&Action', 'S');
Cancel = accelerator(Cancel,'&Action', 'Q');
Reset = accelerator(Reset,'&Action', 'Z');
%==place buttons in panel
Act = add(Act,OK, 'below', [0 -3]);
Act = add(Act,Cancel, nextto(OK), [10 0]);
Act = add(Act,Reset,nextto(Cancel), [10 0]);
Act = marginalize(Act,[3 5]);

%==================
function local_plot(figh);
% close existing plot (if any), launch new plot, and store handle to plot figure
P = getGUIdata(figh,'currentPTC');
plotfigh = getGUIdata(figh,'PlotHandle',[]);
if isSingleHandle(plotfigh), delete(plotfigh); end
ah = plot(P);
plotfigh =  parentfigh(ah(1));
set(plotfigh, 'HandleVisibility', 'off'); % prevent others from accidentily plotting into this figure
setGUIdata(figh,'PlotHandle',plotfigh);
figure(figh); % bring GUI to the fore

%===================
function local_delete(figh); 
% delete the edit figure and the plot, if it exists
plotfigh = getGUIdata(figh,'PlotHandle',[]);
if isSingleHandle(plotfigh), delete(plotfigh); end
delete(figh);

%===================
function local_cmd(figh); 
% get command from GUI
S = GUIval(figh);
if isempty(S), return; end % error has been dealt with by GUIval
% parse command & give feedback in case of trouble
Cmd = Words2cell(S.Cmd);
AllVerbs = {'smooth' 'clip'};
okay=0; % pessimistic default
if numel(Cmd)<3, 
    GUImessage(figh, {'Command must have syntax:' 'f1 f2 Foo param1 ...'}, 'error', 'Cmd');
elseif ~ismember(lower(Cmd{3}),AllVerbs),
    GUImessage(figh, [{['Unknown command ''' Cmd{3} '''.'] 'Known commands are:'} prettyprint(AllVerbs)], 'error', 'Cmd');
elseif isnan(str2double(Cmd{1})),
    GUImessage(figh, 'First word of command must be lower frequency (kHz)', 'error', 'Cmd');
elseif isnan(str2double(Cmd{2})),
    GUImessage(figh, 'Second word of command must be lower frequency (kHz)', 'error', 'Cmd');
else,
    okay=1;
end
if ~okay, return; end
% try to execute the command
P = getGUIdata(figh,'currentPTC');
FreqRange = sort(1e3*[str2double(Cmd{1}) str2double(Cmd{2})]); % freq range in Hz
fun = lower(Cmd{3});
try,
    P = feval(fhandle(fun), P, FreqRange, Cmd{4:end});
catch
    GUImessage(figh, {['Error from fcn' upper(Cmd{3}) ':'], lasterr}, 'error', 'Cmd');
    return;
end
setGUIdata(figh,'currentPTC', P);
local_plot(figh);






