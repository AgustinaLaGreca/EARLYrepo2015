function PlayRecordStop(Mode, hdashboard, hstim, Exp, Recording);
% PlayRecStop - play/record/stop stimulus & data recording
%   PlayRecordStop(hdashboard, hstim, Mode, Exp, RecordingParam, RecParam)
%   launches a Play, PlayRecord, or Stop action. Inputs:
%             Mode: one of 'Play', 'PlayRec', 'Stop'
%       hdashboard: handle to dashboard GUI (needed for messaging on DA progress)
%            hstim: handle to stimulus GUI (needed for stimulus waveforms & info)
%              Exp: experiment to which this action belongs 
%        Recording: struct containing settings on recording returned by
%                   RecordingInstructions, combined with circuit info and 
%                   additional settings from dashboard GUI.
%
%   When Mode='Stop', only the hdashboard arg is needed.
%   At the time of this call, the stimulus must have been checked (see
%   stimGUI). 
%
%   PlayRecordStop is typically called by Dashboard when one of its Action
%   buttons is pushed. It is possible, however, to call PlayRecordStop from
%   the commandline, a function, or script.
%
%   See also Dashboard, RecordingInstructions.

persistent Stim WavPtr % needed for Replay


switch(Mode),
    case {'Play' 'PlayRecord' 'Replay'},
        more off; % when debugging, halted screen output causes problems
        doRec = isequal('PlayRecord', Mode);
        doReplay = isequal('Replay', Mode); 
        if ~doReplay, % first play action: proceed from scratch
            delete(timerfindall); % kill any background action
            Stim = getGUIdata(hstim,'StimParam');
            seqplayinit(Stim.Fsam/1e3, Exp.AudioDevice, '#noload', 1); % '#noload': don't load the circuit; it was already loaded by loadCircuits, called by local_DA below
            WavPtr = GUIaccess(hstim, 'StimParam.Waveform'); % this is where the waveforms are (use a pointer to save memory)
        end
        
        % the first action in the line-up sets the dashboard GUImode to
        % Play or RecPlay when started; put in front because other actions tend to interfere with it.
        actabit(hdashboard, 'setGUImode1', 'start', {@dashboard 'GUImode' Mode});
        % prepare data grabbing and storage (must precede D/A prep because timing calibration uses the DAC hardware)
        if doRec, % prepare time stamping, data transfer and data storage
            % initialize dataset; tell it what the stimulus is
            Pref = preferences(current(experiment));
            if isequal(Pref.CheckFullDsInfo, 'Yes'), 
                CheckFullDsInfo = 1; 
                StE = status(Exp);
                ds.ID.iRecOfCell = StE.iRecOfCell+1;
                ds.ID.iCell = max(1, StE.iCell);
                ds.StimType = Stim.StimType;
                ds_info = PlayRecordSeqID(ds);
                clear ds  StimParam;

            end
            DS = dataset(Exp, Stim, Recording, hdashboard, 'Dataset'); % uploaded
            if exist('CheckFullDsInfo') 
                if ~isempty(ds_info)
                    DS = AddDsInfo(DS,ds_info);
                    upload(DS);
                else
                    local_GUImode(hdashboard, 'Ready');
                    GUImessage(hdashboard, 'Recording Interrupted');
                    return;
                end

            end
            % add the dataset to the dashboard
            setGUIdata(hdashboard,'dataset',DS);
            % launch data grabbers one by one
            for ii=1:numel(Recording.RecordInstr),
                rp = Recording.RecordInstr(ii);
                rp.grabber(rp, rechardware(Exp), DS, hdashboard, ['grab_from_' rp.datafieldname], 'Playit');
            end
            % make sure that data will get completed & saved. DS takes initiative once grabbing is over.
            actabit(hdashboard, 'finalizeDS', 'wrapup', {@getdata DS}, 'clear', {@save DS});
            % online data analysis. Placed behind finalizeDS so that last plot gets complete data
            savedataprogress(DS, 10000, Stim.Dataviewer, hdashboard, 'ViewReponse', 'Playit',Stim.DataviewerParamfile);
            %showdataprogress(DS, 10000, Stim.Dataviewer, hdashboard, 'ViewResponse', 'Playit');
            %showdataprogress(DS, 10e3, @cyclehisto, figh, 'ViewResponse', 'Playit');
        else,% i.e., Play, don't record. This is always repetitive.
            % Original code
            actabit(hdashboard, 'repeatplay', 'clear', {@local_replay Mode, hdashboard, hstim, Exp, Recording});
        end
        % initialize the act of playing the sound
        playsound(WavPtr, Stim.Experiment, Recording, Stim.Attenuation, Stim.Presentation, hdashboard, 'Playit');
        % initialize real-time feedback during D/A
        if exist('DS','var')
            ds_name = IDstring(DS,'status');
        else
            ds_name = '';
        end
        reportstimprogress(Stim.Presentation, hdashboard, 100, hdashboard, 'StimReporter', 'Playit',ds_name,Stim); % playit ready -> stimreporter ready
        % penultimate action in the line-up resets GUImode to busy & ready when the time is there
        actabit(hdashboard, 'setGUImode2', ...
            'wrapup', {@dashboard 'GUImode' 'Busy'}, 'clear', {@dashboard 'GUImode' 'Ready'});
        % the stopper that handles interrupts is last, because it must stop the others before stopping itself.
        stopguiaction(hdashboard, 'Stopper', 'Playit');
        %LL=GUIactionList(hdashboard); {LL.FN}
        % --------------------------------------
        % Let all action objects in hdashboard prepare themselves in order of appearence
        GUIaction(hdashboard, @prepare);
        % Now the action starts
        GUIaction(hdashboard, @start, @wrapup, @clear); % invoke wrapup as soon all action has stopped, then clear
    case 'Stop',
        local_replay('Stop'); % suppress looping Play
        %retrieve Stopper object. On fast 2x clicking Stop button, the stopper may be cleared at the 2nd click, therefore ... 
        Stopper = getGUIdata(hdashboard, 'Stopper',[]); % ...to avoid the error of a vanished Stopper, return [] by default.
        if ~isempty(Stopper), setflag(Stopper);  end % setflag is also made robust against meantime disappearance of Stopper.
end


%=================================
function local_replay(Mode, hdashboard, hstim, Exp, Recording);
persistent StopIt

if isempty(StopIt), StopIt=0; end

if isequal('Stop', Mode),
    StopIt = 1;
    return;
end
if StopIt, % do nothing now, but reset stopIt flag
    StopIt = 0; 
else, % launch repeat play
    ThandlE= timer('TimerFcn', @(Src,Ev)PlayRecordStop('Replay', hdashboard, hstim, Exp, Recording), ...
        'ExecutionMode', 'singleShot');
    start(ThandlE);
end

function local_GUImode(figh, Mode, gmess, gmessmode);
% set enable status of dashboard uicontrols
if nargin<3, gmess = inf; end % indicates absence of message - '' would be bad choice.
if nargin<4, gmessmode = 'neutral'; end
[Mode, Mess] = keywordMatch(Mode,{'Busy' 'Ready' 'Play' 'PlayRecord' 'Replay' 'Stop'}, 'Mode argument');
error(Mess);
A = getGUIdata(figh, 'ActionButton');
A_ExpStat = A('NewUnit', 'NewElectrode', 'Note', 'FinishExp', 'EditExp'); % buttons changing  the Experiment status
Q = getGUIdata(figh, 'Query', ParamQuery()); if numel(Q) == 1 && isvoid(Q), Q = Q([]); end
Exp = current(experiment);
[dum RecSrc] = recordingsources(Exp);
Qmeasure = Q(RecSrc{:}); % queries having to do w recordings
hasExp = ~isvoid(Exp); % true if experiment is going on
if isvoid(Exp), ExpStr = ' (no experiment)'; else, ExpStr = [' -- Experiment ' name(Exp) ' by ' experimenter(Exp)]; end
set(figh, 'name', ['Dashboard' ExpStr]);

hasStim = isSingleHandle(getGUIdata(figh, 'StimGUIhandle'));
switch Mode,
    case {'Busy', 'Stop'}, % disable all buttons & prevent closing the figure
        enable(A,0); enable(Q,0);
        GUIclosable(figh,0); % not okay to close GUI
        % color Check or Stop buttons
        if isequal('Stop', Mode),
            highlight(A('Stop'),[0.5 0.15 0]);
        end
    case 'Ready', % enable all buttons except Stop; okay to close GUI; recording queries depend on experiment status
        enable(A,1);  enable(Q,1); enable(Qmeasure, hasExp);
        %enable(A('StimSpec'), hasExp);
        enable(A('Stop'),0);
        % only enable Play if D/A is possible; only enable PlayRec when an experiment is ongoing
        enable(A('Play'), CanPlayStim && hasStim);
        enable(A('PlayRec'), CanPlayStim && hasStim && canrecord(Exp));
        enable(A_ExpStat, hasExp); % Exp status may only be changed when an Exp has been defined
        enable(A('NewExp', 'ResumeExp'), ~hasExp);
        highlight(A,'default');
        GUIclosable(figh,1); % okay to close GUI
    case {'Play' 'PlayRecord' 'Replay'}, % disable all buttons except Stop
        enable(A,0); enable(Q,0);
        enable(A('Stop'),1);
        GUIclosable(figh,0); % not okay to close GUI
        % color Play or PlayRec buttons
        if isequal('Play', Mode) || isequal('Replay', Mode),
            highlight(A('Play'),[0 0.7 0]);
        elseif isequal('PlayRecord', Mode),
            highlight(A('PlayRec'),[0.85 0 0]);
        end
end
% display GUI message, if any.
if ~isempty(gmess) && ~isequal(inf, gmess),
    GUImessage(figh,gmess,gmessmode);
end
% update the Experiment status info
EM = GUImessenger(figh, 'ExpInfo');
EXP = current(experiment);
reportstatus(EXP, EM);
figure(figh);
drawnow;

