function [Y, Z]=GUI(dum, kw, varargin);
% Probetubecalib/GUI - GUI for performing a probe tube calibration
%   GUI(Probetubecalib(), 'launch') launches the ProbeTubeCalibration GUI,
%   which enables the user to specify the parameters and to measure and
%   save the data of this system calibration. Note that the first arg to
%   GUI is a dummy Probetubecalib object used for Methodizing GUI.
%
%   Probetubecalib/GUI also serves as a callback function for Action
%   buttons, etc, of the actual GUI.
%
%   See also Transfer, Methodizing.

if nargin<2, kw = 'launch'; end

GUIname = 'probetubecalibGUI';
switch lower(kw),
    case 'launch', % ProbeTubeCalibration launch <1>
        Y = local_launch();
        local_GUImode(Y, 'Ready');
    case 'guimode', % GUI(dum, 'guimode', Mode, gmess, gmessmode)
        [ok, figh]=existGUI(GUIname);
        if ~ok, error(['No ' mfilename ' GUI rendered']); end
        local_GUImode(figh, varargin{:});
    case 'keypress',
    case {'check'}, % GUIparams=GUI(dum, 'Check')
        [ok, figh]=existGUI(GUIname);
        if ~ok, error(['No ' mfilename ' GUI rendered']); end
        Y = local_check(figh);
    case {'go', 'stop'}, % GUI(dum, 'Go', 'Left')
        LR = varargin{1};
        figh = gcbf;
        figh = figh.Number;
        if ~isequal('Left', LR), return; end; % ignore right-clicks
        blank(getGUIdata(figh,'Messenger')); % empty messages
        local_DA(figh, kw);
    case {'plot'}, % GUI(dum, 'Plot', 'Left')
        LR = varargin{1};
        figh = gcbf;
        figh = figh.Number;
        if ~isequal('Left', LR), return; end; % ignore right-clicks
        blank(getGUIdata(figh,'Messenger')); % empty messages
        local_plot(figh);
    case {'save'}, % GUI(dum, 'Save', 'Left')
        LR = varargin{1};
        figh = gcbf;
        figh = figh.Number;
        if ~isequal('Left', LR), return; end; % ignore right-clicks
        blank(getGUIdata(figh,'Messenger')); % empty messages
        local_save(figh);
    case {'edit'}, % GUI(dum, 'Edit', 'Left')
        LR = varargin{1};
        figh = gcbf;
        figh = figh.Number;
        if ~isequal('Left', LR), return; end; % ignore right-clicks
        blank(getGUIdata(figh,'Messenger')); % empty messages
        local_edit(figh);
    case {'load'}, % GUI(dum, 'Load', 'Left')
        LR = varargin{1};
        figh = gcbf;
        figh = figh.Number;
        if ~isequal('Left', LR), return; end; % ignore right-clicks
        blank(getGUIdata(figh,'Messenger')); % empty messages
        local_load(figh);
    case {'plot'}, % GUI(dum, 'Plot', ...)
        [ok, figh]=existGUI(GUIname);
        if ~ok, error(['No ' mfilename ' GUI rendered']); end
        local_plot(figh);
    case {'viewfile'}, % GUI(dum, 'viewfile', ...)
        local_viewfile;
    case {'help'}, % GUI(dum, 'Help', ...)
        local_help(varargin{:});
    case 'close',
        if nargin<3, figh=gcbf; else, figh=varargin{1}; end
        figh =figh.Number;
        if ~GUIclosable(figh), return; end % ignore close request if not closable
        % close
        GUIclose(figh);
    otherwise,
        error(['Unknown keyword ''' kw '''.']);
end


%=================================================================
%=================================================================
function figh=local_launch(Exp, launchStim);
% launch ProbeTubeCalibration GUI
dum = probetubecalib(); % dummy Probetubecalib for recursive GUI calls
GUIname = 'probetubecalibGUI';
[EE, figh] = existGUI(GUIname);
if EE, % do not open a 2nd instant of the GUI; just return figh and get out
    figure(figh); return;
end
CLR = 0.75*[1 1 1]+0.1*[0.3 -0.5 0.3];
P_stim = local_stimPanel(CLR); % stimulus parameters
P_rec = local_recPanel(CLR); % recording parameters
P_ax = local_actionPanel(CLR); % play/record/stop & messages
View = pulldownmenu('View', '&View');
View = additem(View,'View Probe tube data from file', @(Src,Ev)GUI(dum, 'viewfile'));
Hlp = pulldownmenu('Help', '&Help');
Hlp = additem(Hlp,'How to calibrate a probe tube', @(Src,Ev)GUI(dum, 'help', 'HowTo'));
%======GUI itself===========
PT=GUIpiece('ProbeTubeCalibration',[],[0 0],[10 4]);
PT = add(PT,P_stim);
PT = add(PT,P_rec, below(P_stim));
PT = add(PT, P_ax, nextto(P_rec));
PT = marginalize(PT,[40 20]);
PT = add(PT, View);
PT = add(PT, Hlp);
% help pulldown menu
% open figure and draw GUI in it
figh = newGUI(GUIname, 'Probe Tube Calibration', {fhandle(mfilename), 'launch'}, 'color', CLR);
draw(figh, PT); 
% empty all edit fields & message field
GUImessage(figh, ' ','neutral');
% closereq & keypress fun
set(figh,'keypressf',{@(Src,Ev)GUI(dum, 'keypress')});
set(figh,'closereq',{@(Src,Ev)GUI(dum, 'close')});
% store StimGUI handle in GUIdata (yet empty)
setGUIdata(figh,'ProbeTubeTrf', []);
% restore most recent, approved settings
GUIfill(figh,0);

%======
function S = local_stimPanel(CLR); 
% panel for specifying stimulus params
S = GUIpanel('Stim', 'stimulus', 'backgroundcolor', CLR);
PeakAmp = ParamQuery('PeakAmp',  'D/A ampl:', '0.125XXX', 'V', 'positive/rreal',  'Peak amplitude of voltage sent to D/A.', 1);
Atten = ParamQuery('Atten',  'Attenuate:', '37.8  ', 'dB', 'nonnegative/rreal',  'Analog attenuation of D/A signal.', 1);
DAchan = ParamQuery('DAchan',  'DA channel:', '', {'1' '2'}, '',  'Click button to switch analog D/A channel of stimulus.', 1);
Fmin = ParamQuery('Fmin',  'Min freq:', '9999999 ', 'Hz', 'positive/rreal',  'Lowest frequency.', 1);
Fmax = ParamQuery('Fmax',  'Max freq:', '9999999 ', 'Hz', 'positive/rreal',  'Highest frequency.', 1);
Speed = ParamQuery('Speed',  'FM speed:', '0.666 ', 'octave/s', 'rreal',  'Number of octaves sewpt per second. Negative numbers reult in downward sweeps.', 1);
S = add(S,PeakAmp);
S = add(S,Atten, below(PeakAmp));
S = add(S,DAchan, below(Atten));
S = add(S,Fmin, nextto(PeakAmp),[20 0]);
S = add(S,Fmax, below(Fmin));
S = add(S,Speed, below(Fmax));
S = marginalize(S,[3 5]);

%======
function S = local_recPanel(CLR); 
% panel for specifying recording params
S = GUIpanel('Rec', 'recording', 'backgroundcolor', CLR);
ProbeSens = ParamQuery('ProbeSens',  'Probe mic:', '0.222XX', 'V/Pa', 'positive/rreal',  'Sensitivity of probe microphone + amplifier (AD channel 1).', 1);
CavitySens = ParamQuery('CavitySens',  'Cavity mic:', '0.999XX', 'V/Pa', 'positive/rreal',  'Sensitivity of cavity microphone + amplifier (AD channel 2).', 1);
S = add(S,ProbeSens);
S = add(S,CavitySens, below(ProbeSens));
S = marginalize(S,[3 5]);

%======
function Act = local_actionPanel(CLR);
% Play/PlayRec/Stop panel
dum = probetubecalib(); % dummy Probetubecalib for recursive GUI calls
Act = GUIpanel('Act', 'action', 'backgroundcolor', CLR);
MessBox = messenger('@MessBox', 'The problem is what you think it is or not ...',5, ... % the '@' in the name indicates that ...
    'fontsize', 12, 'fontweight', 'bold'); % MessBox will be the Main Messenger of he GUI
Go = ActionButton('Go', 'GO', 'XXXXXXXX', 'Start Recording.', @(Src,Ev,LR)GUI(dum, 'Go', LR));
Stop = ActionButton('Stop', 'STOP', 'XXXXXXXX', 'Interrupt ongoing recording.', @(Src,Ev,LR)GUI(dum, 'Stop', LR));
Plot = ActionButton('Plot', 'PLOT', 'XXXXXXXX', 'Plot measured probe-tube transfer.', @(Src,Ev,LR)GUI(dum, 'Plot', LR));
Save = ActionButton('Save', 'SAVE', 'XXXXXXXX', 'Save measured probe-tube transfer.', @(Src,Ev,LR)GUI(dum, 'Save', LR));
Edit = ActionButton('Edit', 'EDIT', 'XXXXXXXX', 'Filter and/or constrain the TRF function.', @(Src,Ev,LR)GUI(dum, 'Edit', LR));
Load = ActionButton('Load', 'LOAD', 'XXXXXXXX', 'Load probe-tube transfer data from file and proceed as if these have just been measured.', @(Src,Ev,LR)GUI(dum, 'Load', LR));
%==fast keys
Go = accelerator(Go,'&Action', 'G');
Stop = accelerator(Stop,'&Action', 'W');
Plot = accelerator(Plot,'&Action', 'P');
Save = accelerator(Save,'&Action', 'S');
Edit = accelerator(Edit,'&Action', 'E');
Load = accelerator(Load,'&Action', 'O');
Act = add(Act,MessBox);
%==place buttons in panel
Act = add(Act,Go, 'below', [0 -3]);
Act = add(Act,Stop,nextto(Go), [10 0]);
Act = add(Act,Plot, nextto(Stop), [10 0]);
Act = add(Act,Save,below(Go), [0 5]);
Act = add(Act,Edit,nextto(Save), [10 0]);
Act = add(Act,Load,nextto(Edit), [10 0]);
Act = marginalize(Act,[3 5]);

%======
function local_GUImode(figh, Mode, gmess, gmessmode);
% set enable status of ProbeTubeCalibration uicontrols
if nargin<3, gmess = inf; end % indicates absence of message - '' would be bad choice.
if nargin<4, gmessmode = 'neutral'; end
[Mode, Mess] = keywordMatch(Mode,{'Busy' 'Ready' 'Play' 'Stop'}, 'Mode argument');
error(Mess);
A = getGUIdata(figh, 'ActionButton');
Q = getGUIdata(figh, 'Query');
TRF_measured = ~isempty(getGUIdata(figh,'ProbeTubeTrf',[]));
switch Mode,
    case {'Busy', 'Stop'}, % disable all buttons & prevent closing the figure
        enable(A,0); enable(Q,0);
        GUIclosable(figh,0); % not okay to close GUI
        % color Check or Stop buttons
        if isequal('Stop', Mode),
            highlight(A('Stop'),[0.5 0.15 0]);
        end
    case 'Ready', % enable all buttons except Stop; okay to close GUI; recording queries depend on experiment status
        enable(A,1);  enable(Q,1); 
        enable(A('Stop'),0);
        % only enable Play if D/A is possible; only enable PlayRec when an experiment is ongoing
        enable(A('Go'), CanPlayStim);
        enable(A('Plot', 'Save', 'Edit'), TRF_measured);
        highlight(A,'default');
        GUIclosable(figh,1); % okay to close GUI
    case {'Play'}, % disable all buttons except Stop
        enable(A,0); enable(Q,0);
        enable(A('Stop'),1);
        GUIclosable(figh,0); % not okay to close GUI
        % color Play or PlayRec buttons
        highlight(A('Go'),[0 0.7 0]);
end
% display GUI message, if any.
if ~isempty(gmess) && ~isequal(inf, gmess),
    GUImessage(figh,gmess,gmessmode);
end
figure(figh);
drawnow;

%======
function S = local_check(figh);
% read stim parameters and check them
setGUIdata(figh, 'ProbeTubeTrf',[]); % erase any previous params
local_GUImode(figh, 'Busy', 'Checking parameters...'); % will be changed to Play or Record inside PlayRecordStop call
S = GUIval(figh); % read parameters. [] is returned when something is wrong.
if ~isempty(S), % params okay - store them in GUI history
    GUIgrab(figh, '>');
end
local_GUImode(figh, 'Ready', ' '); 

%======
function local_DA(figh, kw);
% D/A -related action: Play, PlayRecord or Stop.
% = check stimulus params
dum = probetubecalib(); % dummy Probetubecalib for recursive GUI calls
switch kw,
    case {'Go'}, % prepare D/A
        Param= GUI(dum, 'check');
        if isempty(Param), return; end;
        % perform the calibration
        local_GUImode(figh, 'Play', 'Playing/recording ...'); 
        [Trf, Mess, ErrField] = measure(probetubecalib(),Param, figh);
        GUImessage(figh, Mess, 'error', ErrField);
        if isvoid(Trf), % interrupted
            local_GUImode(figh, 'Ready');
            return; 
        end 
        setGUIdata(figh, 'ProbeTubeTrf', Trf);
        local_GUImode(figh, 'Ready', 'Measurement completed.'); 
        local_plot(figh);
    case 'Stop',
        local_GUImode(figh, 'Stop', 'Stopping ... .', 'warn'); 
        measure(probetubecalib(),'stop');
        seqplayhalt;
        local_GUImode(figh, 'Ready', 'Measurement interrupted.', 'warn'); 
        %PlayRecordStop(kw, figh);
    otherwise,
        error(['Invalid keyword ''' kw '''.']);
end % switch/case

%======
function local_plot(figh);
% plot trf measured functions
figure('numbertitle','off','integerhandle', 'off', 'name', 'Measured Probe Tube Calibration.');
Trf = getGUIdata(figh, 'ProbeTubeTrf');
plot(Trf);
%======
function local_save(figh);
Trf = getGUIdata(figh, 'ProbeTubeTrf');
[okay, FFN] = save(Trf, '?', 'Specify filename for Probe Tube Transfer');
if okay,
    GUImessage(figh, {'Probe Tube calibration data saved to file' FFN });
else,
    GUImessage(figh, 'Cancelled saving data.', 'warn');
end
%======
function local_edit(figh);
Trf = getGUIdata(figh, 'ProbeTubeTrf');
local_GUImode(figh, 'Busy', 'Editing Probe-tube calibration data ...');
[Trf, wasEdited] = edit(Trf);
if wasEdited,
    setGUIdata(figh, 'ProbeTubeTrf', Trf);
    local_GUImode(figh, 'Ready', 'Probe-tube calibration data Edited.');
else, % editing cancelled
    local_GUImode(figh, 'Ready', 'Cancelled editing of probe-tube calibration data.', 'warn');
end

%======
function local_load(figh);
% load probetubecalib data to GUI as if they have just been measured
Trf = getGUIdata(figh, 'ProbeTubeTrf');
if ~isempty(Trf),
    Prompt = {'Loading data from file will overwrite measured data.' 'Continue anyhow?'};
    Answer=questdlg(Prompt,'Warning: data about to be erased.', 'Overwrite', 'Cancel', 'Cancel');
    if isequal('Cancel', Answer), return; end
end
[Trf, FFN] = load(probetubecalib);
if isvoid(Trf),
    GUImessage('Cancelled loading data.', 'warn');
else,
    setGUIdata(figh, 'ProbeTubeTrf', Trf);
    GUI(probetubecalib,'GUImode','Ready', {'Loaded calib data from file', FFN});
end

%======
function local_viewfile;
% view probetubedata from file
[PTC, FFN] = load(probetubecalib());
[dum, Nam] = fileparts(FFN);
if ~isvoid(PTC),
    figure('numbertitle','off','integerhandle', 'off', 'name', [Nam ' (probe tube calibration).']);
    plot(PTC);
end

%======
function local_help(Topic);
switch lower(Topic),
    case 'howto',
        winopen(which('HowToCalibrateAProbeTube.pdf'));
    otherwise,
        error(['Unknown help topic ''' Topic '''.']);
end
    
%=========









