function [Y, Z]=GUI(dum, kw, varargin);
% Earcalib/GUI - GUI for performing an in situ ear calibration
%   EC = GUI(earcalib(), 'launch') launches the earcalib GUI,
%   which enables the user to specify the parameters and to measure and
%   save the data of this in situ calibration. Note that the first input 
%   argument is a dummy Earcalib object serving to Methodize this function.
%
%   Earcalib/GUI also serves as a callback function for Action
%   buttons, etc, of the actual GUI.
%
%   See also Earcalib, Transfer, Methodizing.

if nargin<2, kw = 'launch'; end

GUIname = 'earcalibGUI';
switch lower(kw),
    case 'launch', % earcalibration launch DAchan
        if isvoid(current(experiment)),
            %by abel make this into a warning to resolve dashboard crash
            warning('No current experiment defined. Initialize an experiment before measuring an Earcalib.');
            Y = [];
            Z = [];
            return;
        end
        figh = local_launch(dum);
        local_GUImode(figh, 'Ready');
%         figh = figh.Number;
        GUIreaper(figh, 'wait'); % wait for a trigger to reap ...
        Y = getGUIdata(figh, 'qqqqq',[]); % ... grab the data ...
        local_delete(figh); % ... delete the figure & dependables
    case 'guimode', % GUI(dum, 'guimode', Mode, gmess, gmessmode)
        [ok, figh]=existGUI(GUIname);
        if ~ok, error(['No ' mfilename ' GUI rendered']); end
        local_GUImode(figh, varargin{:});
    case 'keypress',
    case {'filesave' 'fileopen' 'fileemptylist'},
        figh = gcbf;
        figh = figh.Number;
        local_file(figh, kw, varargin);
    case {'probenamebutton'}, % GUI(dum, probenamebutton, LR, chan)
        LR = varargin{1}; 
        if ~isequal('Left', LR), return; end; % ignore right-clicks
        local_PNbutton(gcbf, varargin{2});
    case {'speakernamebutton'}, % GUI(dum, probenamebutton, LR, chan)
        warning('speakernamebutton NYI')
    case {'check'}, % GUIparams=GUI(dum, 'Check')
        [ok, figh]=existGUI(GUIname);
        if ~ok, error(['No ' mfilename ' GUI rendered']); end
        Y = local_check(figh);
    case {'go', 'stop'}, % GUI(dum, 'Go', 'Left', ichan)
        LR = varargin{1};
        figh = gcbf;
        figh = figh.Number;
        if ~isequal('Left', LR), return; end; % ignore right-clicks
        blank(getGUIdata(figh,'Messenger')); % empty messages
        local_DA(figh, kw, varargin{2:end});
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
        GUIreaper(figh, 'reap'); % triggers deletion @ hlated launch
    otherwise,
        error(['Unknown keyword ''' kw '''.']);
end


%=================================================================
%=================================================================
function figh=local_launch(EC);
% launch earcalibration GUI
dum = earcalib(); % dummy earcalib for recursive GUI calls
GUIname = 'earcalibGUI';
[EE, figh] = existGUI(GUIname);
if EE, % do not open a 2nd instant of the GUI; just return figh and get out
    figure(figh); return;
end
CLR = 0.75*[1 1 1]+0.1*[1 0 -1];
P_stim = local_stimPanel(CLR); % stimulus parameters
P_rec = local_recPanel(CLR); % recording parameters
P_ax = local_actionPanel(CLR); % play/record/stop & messages
View = pulldownmenu('View', '&View');
View = additem(View,'View Probe tube data from file', @(Src,Ev)GUI(dum, 'viewfile'));
Hlp = pulldownmenu('Help', '&Help');
Hlp = additem(Hlp,'How to perform an in situ calibration', @(Src,Ev)GUI(dum, 'help', 'HowTo'));
Pfile = local_FileMenu; % File pulldown menu
%======GUI itself===========
PT=GUIpiece('earcalibration',[],[0 0],[10 4]);
PT = add(PT,Pfile);
PT = add(PT,P_rec);
PT = add(PT,P_stim, below(P_rec),[30 0]);
PT = add(PT, P_ax, below(P_stim),[60 0]);
PT = marginalize(PT,[40 20]);
PT = add(PT, View);
PT = add(PT, Hlp);
% open figure and draw GUI in it
figh = newGUI(GUIname, ['Ear Calibration for experiment ' name(current(experiment))], ...
    {fhandle(mfilename), EC, 'launch'}, 'color', CLR);
draw(figh, PT); 
% empty all edit fields & message field
GUImessage(figh, ' ','neutral');
% closereq & keypress fun
set(figh,'keypressf',{@(Src,Ev)GUI(dum, 'keypress')});
set(figh,'closereq',{@(Src,Ev)GUI(dum, 'close')});
% store GUI "globals"
EC.Experiment = current(experiment);
DAchan = EC.Experiment.AudioChannelsUsed;

setGUIdata(figh,'DAchan', DAchan);
setGUIdata(figh,'Earcalib', EC);
setGUIdata(figh,'Trf1', transfer()); % future transfer of AD chan 1
setGUIdata(figh,'Trf2', transfer()); % future transfer of AD chan 2
% restore most recent, approved settings
GUIfill(figh,0);
local_set_probes(EC,figh);


function local_set_probes(EC,figh)
figure(figh);
Probes = EC.Experiment.Audio.Probes;
Speakers = EC.Experiment.Audio.Speakers;

Q = getGUIdata(figh,'Query',[]);
GSQ = GUIsettings('ParamQuery');
for i=1:numel(Q)
    q = Q(i);
    if strcmpi(name(q), 'ProbeName_1')
        EditStr = Probes{1};
        he = getFieldOrDefault(q.uiHandles, 'Edit', []);
        if ~isempty(EditStr) && isUIcontrol(he), % impose string value on q 
            set(he,'string', EditStr, 'backgroundcolor', GSQ.EditOkayColor);
        end
    end
    if strcmpi(name(q), 'SpeakerName_1')
        EditStr = Speakers;
        he = getFieldOrDefault(q.uiHandles, 'Edit', []);
        if ~isempty(EditStr) && isUIcontrol(he), % impose string value on q 
            set(he,'string', EditStr, 'backgroundcolor', GSQ.EditOkayColor);
        end
    end
    if strcmpi(name(q), 'ProbeName_2')
        EditStr = Probes{2};
        he = getFieldOrDefault(q.uiHandles, 'Edit', []);
        if ~isempty(EditStr) && isUIcontrol(he), % impose string value on q 
            set(he,'string', EditStr, 'backgroundcolor', GSQ.EditOkayColor);
        end
    end 
    if strcmpi(name(q), 'SpeakerName_2')
        EditStr = Speakers;
        he = getFieldOrDefault(q.uiHandles, 'Edit', []);
        if ~isempty(EditStr) && isUIcontrol(he), % impose string value on q 
            set(he,'string', EditStr, 'backgroundcolor', GSQ.EditOkayColor);
        end
    end
end

%======
function Pfile = local_FileMenu; 
% File pulldown menu
dum = earcalib(); % dummy earcalib for recursive GUI calls
Pfile=pulldownmenu('File','&File');
Pfile=additem(Pfile,'&Save parameters', @(Src,Ev)GUI(dum, 'FileSave',Src), 'accelerator', 'S');
Pfile=additem(Pfile,'&Open parameter file', @(Src,Ev)GUI(dum, 'FileOpen',Src), 'accelerator', 'O');
Pfile=additem(Pfile,'&Empty File List', @(Src,Ev)GUI(dum, 'FileEmptylist',Src));
CL=cycleList('readagain', 10, @(Src, Event, L_Callback, L, Item)GUI(dum,'FileReadagain',L_Callback, L, Item)); % clicking an item calls this fcn with 4 input args; see CycleList
Pfile=additem(Pfile,CL);
%======
function S = local_stimPanel(CLR); 
% panel for specifying stimulus params
S = GUIpanel('Stim', 'stimulus', 'backgroundcolor', CLR);
Fmin = ParamQuery('Fmin',  'Min freq:', '9999999 ', 'Hz', 'positive/rreal',  'Lowest frequency.', 1);
Fmax = ParamQuery('Fmax',  'Max freq:', '9999999 ', 'Hz', 'positive/rreal',  'Highest frequency.', 1);
Speed = ParamQuery('Speed',  'FM speed:', '0.666 ', 'octave/s', 'rreal',  'Number of octaves swept per second. Negative numbers result in downward sweeps.', 1);
MaxSPL = ParamQuery('MaxSPL',  'Max SPL:', 'XXXX', 'dB SPL', 'rreal',  'Maximum SPL at eardrum during calibration.', 1);
DAchan = ParamQuery('DAchan',  'DA channel:', '', {'1' '2'}, '',  'Click button to switch analog D/A channel of stimulus.', 1);
S = add(S,Fmin);
S = add(S,Fmax, nextto(Fmin), [15 0]);
S = add(S,Speed, nextto(Fmax), [15 0]);
S = add(S,MaxSPL, nextto(Speed), [20 0]);
S = marginalize(S,[3 5]);

%======
function S = local_recPanel(CLR); 
% panel for specifying recording params
dum = earcalib(); % dummy earcalib for recursive GUI calls
S = GUIpanel('Rec', 'Probes & Speakers', 'backgroundcolor', CLR);
ProbeName_1 = ParamQuery('ProbeName_1',  'Probe tube 1:', 'XXXXXXXXXXXXXXXXXXXXXXX', '', 'string',  'Name of Probe Tube @ AD channel 1.', 1e2);
ProbeName_2 = ParamQuery('ProbeName_2',  'Probe tube 2:', 'XXXXXXXXXXXXXXXXXXXXXXX', '', 'string',  'Name of Probe Tube @ AD channel 2.', 1e2);
ProbeNameButton_1 = ActionButton('ProbeNameButton_1',  'Browse...', 'xxxxxxxx', 'Click to browse for a probe tube file.', @(Src,Ev,LR)GUI(dum, 'ProbeNameButton', LR, 1));
ProbeNameButton_2 = ActionButton('ProbeNameButton_2',  'Browse...', 'xxxxxxxx', 'Click to browse for a probe tube file.', @(Src,Ev,LR)GUI(dum, 'ProbeNameButton', LR, 2));
SpeakerName_1 = ParamQuery('SpeakerName_1',  'Speaker 1:', 'XXXXXXXXXXXXXXX', '', 'string',  'Name of loudspeaker @ AD channel 1.', 1e2);
SpeakerName_2 = ParamQuery('SpeakerName_2',  'Speaker 2:', 'XXXXXXXXXXXXXXX', '', 'string',  'Name of loudspeaker @ AD channel 2.', 1e2);
SpeakerNameButton_1 = ActionButton('SpeakerNameButton_1',  'Browse...', 'xxxxxxxx', 'Click to browse for a loudspeaker.', @(Src,Ev,LR)GUI(dum, 'SpeakerNameButton', LR, 1));
SpeakerNameButton_2 = ActionButton('SpeakerNameButton_2',  'Browse...', 'xxxxxxxx', 'Click to browse for a loudspeaker.', @(Src,Ev,LR)GUI(dum, 'SpeakerNameButton', LR, 2));
ProbeSens_1 = ParamQuery('ProbeSens_1',  'Probe mic:', '0.222XX', 'V/Pa', 'positive/rreal',  'Current setting (!) of sensitivity of probe microphone + amplifier (AD channel 1).', 1);
ProbeSens_2 = ParamQuery('ProbeSens_2',  'Probe mic:', '0.222XX', 'V/Pa', 'positive/rreal',  'Current setting (!) of sensitivity of probe microphone + amplifier (AD channel 2).', 1);
%==

S = add(S, SpeakerName_1,'nextto',[20 0]);
S = add(S, SpeakerNameButton_1, nextto(SpeakerName_1), [-9 0]);
S = add(S, SpeakerName_2, nextto(SpeakerName_1), [170 0]);
S = add(S, SpeakerNameButton_2, nextto(SpeakerName_2), [-9 0]);
S = add(S, ProbeName_1, below(SpeakerName_1),[-15 0]);
S = add(S, ProbeNameButton_1, nextto(ProbeName_1), [-9 0]);
S = add(S,ProbeSens_1, below(ProbeName_1), [15 0]);
S = add(S, ProbeName_2, below(SpeakerName_2),[-15 0]);
S = add(S, ProbeNameButton_2, nextto(ProbeName_2), [-9 0]);
S = add(S,ProbeSens_2, below(ProbeName_2), [15 0]);
S = marginalize(S,[3 5]);

%======
function Act = local_actionPanel(CLR);
% Play/PlayRec/Stop panel
dum = earcalib(); % dummy earcalib for recursive GUI calls
Act = GUIpanel('Act', 'action', 'backgroundcolor', CLR);
MessBox = messenger('@MessBox', 'The problem is what you think it is or not ...',5, ... % the '@' in the name indicates that ...
    'fontsize', 12, 'fontweight', 'bold'); % MessBox will be the Main Messenger of he GUI
Go_1 = ActionButton('Go_1', 'REC 1', 'XXXXXXXX', 'Start Recording over D/A channel 1.', @(Src,Ev,LR)GUI(dum, 'Go', LR, 1));
Go_2 = ActionButton('Go_2', 'REC 2', 'XXXXXXXX', 'Start Recording over D/A channel 2.', @(Src,Ev,LR)GUI(dum, 'Go', LR, 2));
Save = ActionButton('Save', 'SAVE', 'XXXXXXXX', 'Save measured ear calibration.', @(Src,Ev,LR)GUI(dum, 'Save', LR));
Stop = ActionButton('Stop', 'STOP', 'XXXXXXXX', 'Interrupt ongoing recording.', @(Src,Ev,LR)GUI(dum, 'Stop', LR));
Plot = ActionButton('Plot', 'PLOT', 'XXXXXXXX', 'Plot measured ear calibration.', @(Src,Ev,LR)GUI(dum, 'Plot', LR));
Edit = ActionButton('Edit', 'EDIT', 'XXXXXXXX', 'Filter and/or constrain the TRF function.', @(Src,Ev,LR)GUI(dum, 'Edit', LR));
Left = ParamQuery('Left',  'Left', '', '', 'string',  'Left Ear', 1e2);
Right = ParamQuery('Right',  'Right', '', '', 'string',  'Right Ear', 1e2);
%==fast keys
Stop = accelerator(Stop,'&Action', 'W');
Plot = accelerator(Plot,'&Action', 'P');
Save = accelerator(Save,'&Action', 'S');
% Edit = accelerator(Edit,'&Action', 'E');
% Load = accelerator(Load,'&Action', 'O');
Act = add(Act,MessBox,'nextto',[200 0]);
%==place buttons in panel
Act = add(Act,Go_1, below(MessBox),[-200 0]);
Act = add(Act,Left, below(Go_1),[10 0]);
Act = add(Act,Stop,nextto(Go_1), [30 40]);
Act = add(Act,Plot, nextto(Stop), [10 0]);
Act = add(Act,Save,nextto(Plot), [10 0]);
Act = add(Act,Edit,nextto(Save), [10 0]);
Act = add(Act,Go_2, nextto(Edit), [30 -40]);
Act = add(Act,Right, below(Go_2),[10 0]);
Act = marginalize(Act,[3 5]);

%======
function local_GUImode(figh, Mode, gmess, gmessmode);
% set enable status of earcalibration uicontrols
if nargin<3, gmess = inf; end % indicates absence of message - '' would be bad choice.
if nargin<4, gmessmode = 'neutral'; end
[Mode, Mess] = keywordMatch(Mode,{'Busy' 'Ready' 'Play_1' 'Play_2' 'Stop'}, 'Mode argument');
error(Mess);
A = getGUIdata(figh, 'ActionButton');
Q = getGUIdata(figh, 'Query');
% D/A chan specific buttons & queries
ADchan = getGUIdata(figh,'DAchan');
DA1_active = ismember(ADchan,{'Left' 'Both'});
DA2_active = ismember(ADchan,{'Right' 'Both'});
Q_1 = Q('ProbeName_1', 'ProbeSens_1', 'SpeakerName_1'); 
Q_2 = Q('ProbeName_2', 'ProbeSens_2', 'SpeakerName_2'); 
A_1= A('Go_1', 'ProbeNameButton_1', 'SpeakerNameButton_1'); 
A_2= A('Go_2', 'ProbeNameButton_2', 'SpeakerNameButton_2'); 
TRF_measured = ~isempty(getGUIdata(figh,'Earcalib',[]));
EC = getGUIdata(figh, 'Earcalib');
got_Trf1 = ~isvoid(getGUIdata(figh,'Trf1'));
got_Trf2 = ~isvoid(getGUIdata(figh,'Trf2'));
Trfcomplete = (got_Trf1 || ~DA1_active) && (got_Trf2 || ~DA2_active);
Trfabsent = ~got_Trf1 && ~got_Trf2;
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
        % enable buttons & edits according to D/A channels' active states
        enable(A_1, CanPlayStim && DA1_active);
        enable(A_2, CanPlayStim && DA2_active);
        enable(Q_1, CanPlayStim && DA1_active);
        enable(Q_2, CanPlayStim && DA2_active);
        enable(A('Save', 'Edit'),Trfcomplete);
        enable(A('Plot'),~Trfabsent);
        highlight(A,'default');
        if got_Trf1, highlight(A(['Go_1']),[0 0.7 0]); end
        if got_Trf2, highlight(A(['Go_2']),[0 0.7 0]); end
        GUIclosable(figh,1); % okay to close GUI
    case {'Play_1' 'Play_2'}, % disable all buttons except Stop
        enable(A,0); enable(Q,0);
        enable(A('Stop'),1);
        GUIclosable(figh,0); % not okay to close GUI
        % color Play or PlayRec buttons
        if isequal('Play_1', Mode),
            highlight(A(['Go_1']),[0.7 0 0]);
        else,
            highlight(A(['Go_2']),[0.7 0 0]);
        end
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
local_GUImode(figh, 'Busy', 'Checking parameters...');
S = GUIval(figh); % read parameters. [] is returned when something is wrong.
if ~isempty(S), % params okay - store them in GUI history
    GUIgrab(figh, '>');
end
local_GUImode(figh, 'Ready', ' '); 

%======
function local_DA(figh, kw, ichan);
% D/A -related action: Play, PlayRecord or Stop.
% = check stimulus params
dum = earcalib(); % dummy earcalib for recursive GUI calls
switch kw,
    case {'Go'}, % prepare D/A
        chanstr = num2str(ichan);
        Param= GUI(dum, 'check');
        if isempty(Param), return; end;
        % perform the calibration
        local_GUImode(figh, ['Play_' chanstr], 'Playing/recording ...'); 
        EC = getGUIdata(figh, 'Earcalib');
        [Trf, Mess, ErrField] = measure(EC, ichan, ichan, Param, figh);
        if isvoid(Trf), % interrupted or param error
            local_GUImode(figh, 'Ready', Mess, 'warn'); 
            return; 
        end 
        setGUIdata(figh, ['Trf' chanstr], Trf);
        setGUIdata(figh, ['PTC_' chanstr], load(probetubecalib,Param.(['ProbeName_' chanstr])));
        setGUIdata(figh, ['Param' chanstr], Param);
        local_GUImode(figh, 'Ready', 'Measurement completed.'); 
        local_plot(figh);
    case 'Stop',
        local_GUImode(figh, 'Stop', 'Stopping ... .', 'warn'); 
        measure(earcalib(),'stop');
        seqplayhalt;
        local_GUImode(figh, 'Ready', 'Measurement interrupted.', 'warn'); 
        %PlayRecordStop(kw, figh);
    otherwise,
        error(['Invalid keyword ''' kw '''.']);
end % switch/case

%======
function local_plot(figh);
% plot measured TRFs
EC = local_assembleEC(figh);
% attempt to recycle plot figure, but open new one if necessary
hplotfig = getGUIdata(figh,'PlotHandle',[]);	
if ~isSingleHandle(hplotfig),
    hplotfig = figure('integerhandle','off', 'numbertitle', 'off', 'name', 'Measured ear calibration');
    setGUIdata(figh,'PlotHandle',hplotfig);
end
set(hplotfig, 'handlevisib','on');
setGUIdata(figh, 'hplotfig', hplotfig);
figure(hplotfig);
plot(EC);
set(hplotfig, 'handlevisib','off');
figure(figh); % bring GUI to the fore

%======
function EC = local_assembleEC(figh);
% grab all measured data and specified params from GUI and assemble into a
% earcalib object.
EC = getGUIdata(figh, 'Earcalib');
EC.activeDAchan = getGUIdata(figh, 'DAchan');
EC.Param{1} = getGUIdata(figh, 'Param1',[]);
EC.Param{2} = getGUIdata(figh, 'Param2',[]);
PTC_1 = getGUIdata(figh,'PTC_1', []);
PTC_2 = getGUIdata(figh,'PTC_2', []);
% get the probe tube calib data and incorporate them
Trf1 = getGUIdata(figh, 'Trf1',1); Trf2 = getGUIdata(figh, 'Trf2',1);
if ~isvoid(Trf1), 
    Trf1 = Trf1*tubeTrf(PTC_1); 
    Trf1=description(Trf1, '@D/A --> eardrum');
    Trf1 = setWBdelay(Trf1,'wflat');
end
if ~isvoid(Trf2), 
    Trf2 = Trf2*tubeTrf(PTC_2); 
    Trf2=description(Trf2, 'D/A --> eardrum');
    Trf2 = setWBdelay(Trf2,'wflat');
end;
EC.Transfer = [Trf1 Trf2];
% EC.Transfer(1) = Trf1; EC.Transfer(2) = Trf2;

%======
function local_save(figh);
EC = local_assembleEC(figh);
FFN = save(EC);
% remove Trf1 and Trf2, etc, from GUI to prevent saving them twice
setGUIdata(figh, 'Trf1',transfer()); setGUIdata(figh, 'Trf2',transfer());
setGUIdata(figh, 'PTC_1',[]); setGUIdata(figh, 'PTC_2',[]);
FFN = strrep(FFN, datadir, '...');
local_GUImode(figh, 'Ready', {'Ear calib data saved to file' FFN });


%======
function local_edit(figh);
Trf = getGUIdata(figh, 'ProbeTubeTrf');
local_GUImode(figh, 'Busy', 'Editing Probe-tube calibration data ...');
[Trf, wasEdited] = edit(Trf);
if wasEdited,
    setGUIdata(figh, 'ProbeTubeTrf', Trf);
    local_GUImode(figh, 'Ready', 'Probe-tube calibration data Edited.');
else, % editing cancelled
    local_GUImode(figh, 'Ready', 'Cancled editing of probe-tube calibration data.', 'warn');
end

%======
function local_load(figh);
% load earcalib data to GUI as if they have just been measured
Trf = getGUIdata(figh, 'ProbeTubeTrf');
if ~isempty(Trf),
    Prompt = {'Loading data from file will overwrite measured data.' 'Continue anyhow?'};
    Answer=questdlg(Prompt,'Warning: data about to be erased.', 'Overwrite', 'Cancel', 'Cancel');
    if isequal('Cancel', Answer), return; end
end
[Trf, FFN] = load(earcalib);
if isvoid(Trf),
    GUImessage('Cancelled loading data.', 'warn');
else,
    setGUIdata(figh, 'ProbeTubeTrf', Trf);
    GUI(earcalib,'GUImode','Ready', {'Loaded calib data from file', FFN});
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
        winopen(which('HowToCalibrateInSitu.pdf'));
    otherwise,
        error(['Unknown help topic ''' Topic '''.']);
end
    
%=========
function local_PNbutton(figh, ichan);
% browse for probe
chanstr = num2str(ichan);
[P, FFN] = load(probetubecalib(),'?',['Identify probe tube of D/A chan ' chanstr]);
if isvoid(P), return; end
% fill in name in corresponding edit field
figh = figh.Number;
Q = getGUIdata(figh, 'Query');
[dum, FN] = fileparts(FFN);
he = edithandle(Q(['ProbeName_' chanstr]));
set(he, 'string', FN); drawnow;

%===================
function local_delete(figh); 
% delete the plot, if it exists, then the GUI figure
plotfigh = getGUIdata(figh,'PlotHandle',[]);
if isSingleHandle(plotfigh), delete(plotfigh); end
delete(figh);

%===================
function local_file(figh, kw, varargin);
%figh, kw, varargin{:};
switch lower(kw),
    case 'fileopen',
        [Mess, Ignored, FN]=GUIfill(figh,'?',1); % 1: do display message
        [dum FN] = fileparts(FN);
    case 'filesave',
        S=GUIgrab(figh, '?',1);% 1: do display message
    case 'fileemptylist',
        warning(['fileemptylist NIY'])
end % switch/case
local_GUImode(figh,'Ready');






