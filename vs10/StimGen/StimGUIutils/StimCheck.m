function okay = StimCheck(figh, varargin);
% StimCheck - check stimulus parameters of stimulus GUI  
%    StimCheck(figh), where figh is the handle of the stimulus GUI,
%    reads the stimulus parameter from the stimulus GUI and checks their
%    validity. StimCheck returns True if all is well, False otherwise. In
%    the latter case, feedback is provided to the stimulus GUI through its
%    messengers.
%    Next, StimCheck calls the stimulus generator makestimFOO. This may
%    again reveal invalid parameter values or inconsistencies between 
%    parameters. StimCheck will then halt and return False.
%    Finally, StimCheck obtains the generic stimulus parameters by calling
%    GenericStimParams.
%
%    StimCheck is called when the CHECK action button of the stimulus 
%    GUI is pushed, but also implicitly, after a request to play or plot
%    the stimuli.
%
%    See also GenericStimParams.

% reset all messengers
if nargin < 2  || ~isSingleHandle(figh)
    parent_figh = [];
    kw = '';
elseif nargin >= 2
    parent_figh=varargin{1};
    kw = varargin{2};
end
figh = parentfigh(figh);
figh = figh.Number;
allM = getGUIdata(figh,'Messenger');
reset(allM);
% triggered by Check button but also by Play & PlayRec in dashboard
okay=0; % default in case of premature return
% empty stimparam in userdata
setGUIdata(figh,'StimParam', []);
GUImessage(figh, 'Checking ...', 'neutral');
% read params from GUI; only generic checking of user-specified values
P=GUIval(figh); 
if isempty(P), return; end
P.StimType = getGUIdata(figh, 'StimulusType');
P.Experiment = getGUIdata(figh, 'Experiment');
% Pass params to the GUI's stimulus generator, who will do further checking ... 
StimMaker = getGUIdata(figh,'StimMaker'); %... compute waveforms and attach them to P
P=feval(StimMaker, P);
%Check if stutter is on if so repeat the first condition
P = local_AddStutterInfo(P,parent_figh,kw);
if isempty(P), return; end
P.Waveform = defragment(P.Waveform); % defragment the waveforms within P
P.GenericStimParams = GenericStimparams(P); % get generic stimulus parameters (useful for standard analyses)
error(TestNsam(P.Waveform)); % equal sample counts in DAC channels
% Collect dataviewer info into cells
[P, okay] = local_collectDataview(figh, P);
if ~okay, return; end
% Load dataviewparam
[DVparam, DVPokay] = local_loadDataviewParam(figh, P);
if ~DVPokay, return; end
P.DataviewerParamfile = DVparam;
% If we get here, we're ready to play and/or record. Prepare that.
ReportFsam(figh, P.Fsam/1e3);
ReportRMS(figh, P);
GUImessage(figh,'Parameters okay');
% store stimparam and dataviewparam in userdata
setGUIdata(figh,'StimParam', P);
setGUIdata(figh,'DataviewParam', DVparam);
% store GUI state to history
GUIgrab(figh,'>'); % add to history
GUImessage(figh, 'Parameters okay', 'neutral');
refresh(figh);
okay=1;


%=================================================
function [P2, okay] = local_collectDataview(figh,P)
% Assemble all Dataview related Fields into cells
P2 = P;
okay = 0; % pessimistic default
% first dataviewers...
if isequal('-', P.Dataviewer), % no viewer active, so use default dotraster
%     okay = 1;
%     return;
    P.Dataviewer = {'dotraster'}; 
    P.DataGrazeActive = 0;
    P.DataviewerParamfile = 'def';
elseif isequal('active', P.Dataviewer) % viewer(s) active but not specified
    GUImessage(figh, 'No dataviewer specified.', 'error', {'Dataviewer'});
    return;
else
    P.Dataviewer = regexp(erase(P.Dataviewer, '- '), ' ', 'split');
    if numel(P.Dataviewer) > 2,
        GUImessage(figh, 'Too many dataviewers possibly make online analysis slow.', 'warning', P.handle.Dataviewer);
    end
    P.DataGrazeActive = 1;
end
%... and then their corresponding param files
if isequal('def',P.DataviewerParamfile);
    Pfile = cell(size(P.Dataviewer));
    Pfile(:) = {'def'};
    P.DataviewerParamfile = Pfile;
else
    P.DataviewerParamfile = regexp(P.DataviewerParamfile,' ','split');
end
if ~(size(P.Dataviewer)==size(P.DataviewerParamfile))
    GUImessage(figh,...
        'Dataviewer parameters do not correspond to the specified dataviewers.',...
        P.handle.DataviewerParamfile);
    return; 
end;
P2 = P;
okay = 1;

function [DVP, okay] = local_loadDataviewParam(figh, P)
% try to load dataviewparam from file; return void DVP 
okay=1; % optimistic default
% if isequal('-', P.Dataviewer), % this can't be the case anymore
%     DVP = {dataviewparam()}; % no viewer active; return void dataviewparam object
%     return;
% end
% We now know that all dataview related field are cells (and have the same
% size).
DVP = cell(size(P.Dataviewer)); % output is going to be a struct of dataviewparams
for i=1:numel(P.Dataviewer)
    Viewer = P.Dataviewer{i};
    Pfile = P.DataviewerParamfile{i};
    if isequal('def', lower(Pfile))
        DVP{i} = dataviewparam(Viewer); 
    else  % getting serious: try to load the params from file
        StimType = getGUIdata(figh, 'StimulusType');
        subDir = [StimType '_' Viewer]; % dir to load & save dataviewparam
        Pfile = strtok(Pfile,'.');
        DVP{i} = load(dataviewparam(), Viewer, Pfile, subDir);
        okay = ~isvoid(DVP{i});
        if ~okay
            GUImessage(figh, ...
            {[Viewer ' parameter file ''' subDir '\' P.DataviewerParamfile ''''], ...
            'not found.'}, 'error', P.handle.DataviewerParamfile);
            return;
        end
    end
end

% Stutter repeats the first condition twice
% Implemented by Adriaan Lambrechts in November 2016
function P = local_AddStutterInfo(P,figh,kw)
if isempty(figh) || strcmpi(kw,'play') % Stutter should only be used with Play/Record
    return;
end
Query = GUIval(figh);
StutterOn = Query.Stutter;
if strcmpi(StutterOn,'off')
    return;
else
    % get the current number of conditions
    %Ncond = P.Presentation.Ncond;
    % find all variables in the P struct with Ncond rows
    % make for each of these variables a copy of the first row and insert
    % it at the beginning so to create the stutter effect
%     fieldNames = fieldnames(P);
%     for i=1:length(fieldNames)
%         % We only want numeric variables
%         fn = (fieldNames{i});
%         if ~ischar(P.(fn)) && ~strcmpi(fn,'genericparamscall')
%             if length(P.(fn))== Ncond % In this case repeat first condition
%                 P.(fn) = [P.(fn)(1,:); P.(fn)];
%             end
%         end
%     end
end

% Copy the first wave & Attenuation
% P.Waveform = [P.Waveform; P.Waveform(1,:)];
% P.Attenuation.NumScale = [P.Attenuation.NumScale; P.Attenuation.NumScale(1,:)];
% P.Attenuation.NumGain_dB = [P.Attenuation.NumGain_dB; P.Attenuation.NumGain_dB(1,:)];

% Update the Presentation Struct
P.Presentation = addStutterToPres(P.Presentation);

% Add a Stutter field to the P struct
P.Stutter = 'On';
return


% elseif isequal('def', lower(P.DataviewerParamfile)),
%     DVP = dataviewparam(P.Dataviewer); % default for the current dataviewer 
% else, % getting serious: try to load the params from file
%     StimType = getGUIdata(figh, 'StimulusType');
%     subDir = [StimType '_' P.Dataviewer]; % dir to load & save dataviewparam
%     Pfile = strtok(P.DataviewerParamfile,'.');
%     DVP = load(dataviewparam(), P.Dataviewer, Pfile, subDir);
%     okay = ~isvoid(DVP);
%     if ~okay, 
%         GUImessage(figh, ...
%             {['Dataview parameter file ''' subDir '\' P.DataviewerParamfile ''''], ...
%             'not found.'}, 'error', P.handle.DataviewerParamfile);
%     end
% end  
