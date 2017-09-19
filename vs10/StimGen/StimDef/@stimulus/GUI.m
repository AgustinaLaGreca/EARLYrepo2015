function [figh, Mess]=StimGUI(ST, EXP, varargin);
% Stimulus/GUI - open stimulus GUI
%    StimGUI(ST) opens stimulus GUI for stimulus object ST and
%    returns a handle to it. By convention, the stimulus-specific work is 
%    delegated to the overloaded method stimulus/stimdef. The Uhr example
%    is stimFS/stimdef.
%
%    StimGUI('Foo', EXP) applies experiment definition EXP instead of the 
%    current experiment.
%
%    [h, Mess]=StimGUI(..) suppresses error message due to non-existing
%    stimulus paradigms and returns the error in char Mess instead. 
%
%    StimGUI is also the callback function of the "Check" button of the
%    stimulus GUI.
%
%    See also StimDefDir, stimdefFS, Experiment/edit.

if nargin<2,
    EXP = current(experiment);
    if isvoid(EXP),
        EXP = find(experiment,'StimOnly');
    end
end
if isvoid(EXP),
    error('No valid experiment defined; no default (StimOnly) experiment found, either.');
end


% create param part of GUI. Delegate to paramGUI method of ST object
Params = paramGUI(ST, EXP); % Params = stimdefXXX(EXP);

Pfile = StimGuiFilePulldown;
Efile = StimGuiEditPulldown;


%-----------the whole GUI-----------------
SG=GUIpiece([StimName 'Menu'],[],[0 0],[10 4]);
SG = add(SG, Pfile);
SG = add(SG, Efile);
SG = add(SG, Params);
SG = add(SG, local_Action, 'below', [20 10]);
SG = add(SG, local_Dataview, 'nextto', [10 0]);

SG=marginalize(SG,[0 20]);

% open figure and draw GUI in it
figh = newGUI([upper(ST.StimType) 'menu'], [upper(ST.StimType) ' menu   ' name(EXP)], {fhandle(mfilename), ST, EXP});
setGUIdata(figh, 'Stimulus', ST);
CreateQueryContextMenu(figh,ParamQuery);
draw(figh, SG); 
setGUIdata(figh, 'Experiment', EXP);
% empty all edit fields & message field
clearedit(getGUIdata(figh, 'Query'));
GUImessage(figh, ' ','neutral');
% reload any previous pameters
GUIfill(figh,0,0);

%===================================
function Act = local_Action;
% Check/Play/PlayRec dashboard
Act = GUIpanel('Act', '', 'backgroundcolor', 0.75+[0 0.1 0.05]);
MessBox = messenger('@MessBox', 'The problem is what you think the problem is         ',7, ... % the '@' in the name indicates that ...
    'fontsize', 12, 'fontweight', 'bold'); % MessBox will be the Main Messenger of he GUI
FsamDisp = messenger('FsamDisp', 'Fsam = ***.* kHz', 1, 'FontWeight', 'bold', 'ForegroundColor', [0.25 0.15 0.1]);
Check = ActionButton('Check', 'CHECK', 'XXXXXXXXXX', 'Check stimulus parameters and update information.', @(Src,Ev,LR)StimGUI([],'Check',LR));
Check = accelerator(Check,'&Action', 'K');
Act = add(Act,MessBox);
Act = add(Act,FsamDisp, 'nextto', [32 3]);
Act = add(Act,Check, 'below', [0 88]);
Act = marginalize(Act,[3 5]);

function D = local_Dataview;
% Dataview panel
D = GUIpanel('Dataview', 'data view');
DVlist = ['-' listdataviewer(dataset())]; % list of dataviewers plus '-' entry == none
Viewer = ParamQuery('Dataviewer', 'viewer:', '', DVlist, ...
    '', 'Select viewer for online data analysis by clicking the toggle.');
Pfile = ParamQuery('DataviewerParamfile', 'file:', 'XXXXXXXXXX', '', ...
    'varname', 'Provide filename for online dataviewer parameters. "Def" means: use default parameters.', 128);
EditPV = ActionButton('EditDVparams', 'Edit', '----', 'Edit dataviewer parameter values.', ...
    @(Src,Ev,LR)local_edit(Src, name(Pfile), name(Viewer)));
%Check = accelerator(Check,'&Action', 'K');
D = add(D,Viewer);
D = add(D,Pfile,'below');
D = add(D,EditPV,'nextto',[-5 0]);
D = marginalize(D,[0 5]);


function Mess=local_CheckStimfile(fn, StimName, Descr);
% check whether mfile exists, is unique, and has correct location
Mess = '';
FN = fullfile(stimDefDir, [fn '.m']);
if ~exist(FN,'file'),
    Mess = [Descr ' of ''' StimName ''' stimulus not found. See StimGUI for naming conventions.'];
end
qq = which(fn,'-all');
if numel(qq)>1,
    Mess = ['Multiple ' Descr 's present for ''' StimName ''' stimulus. See StimGUI for naming conventions.'];
end

function okay = local_check(figh, varargin);
% ====called when the CHECK baction button of the stimulus GUI is pushed.
% reset all messengers
allM = getGUIdata(figh,'Messenger');
reset(allM);
% triggered by Check button but also by Play & PlayRec in dashboard
okay=0; % default in case of premature return
% empty stimparam in userdata
setGUIdata(figh,'StimParam', []);
GUImessage(figh, 'Checking ...', 'neutral');
% read params from GUI; only generic checking of user-specified values
P=GUIval(figh); if isempty(P), return; end
P.StimType = getGUIdata(figh, 'StimulusType');
P.Experiment = getGUIdata(figh, 'Experiment');
% Pass params to the GUI's stimulus generator, who will do further checking ... 
StimMaker = getGUIdata(figh,'StimMaker'); %... compute waveforms and attach them to P
P=feval(StimMaker, P);
if isempty(P), return; end
P.Waveform = defragment(P.Waveform); % defragment the waveforms within P
error(TestNsam(P.Waveform));
% Load dataviewparam
[DVparam, DVPokay] = local_loadDataviewParam(figh, P);
if ~DVPokay, return; end
% If we get here, we're ready to play and/or record. Prepare that.
ReportFsam(figh, P.Fsam/1e3);
GUImessage(figh,'Parameters okay');
% store stimparam and dataviewparam in userdata
setGUIdata(figh,'StimParam', P);
setGUIdata(figh,'DataviewParam', DVparam);
% store GUI state to history
GUIgrab(figh,'>'); % add to history
GUImessage(figh, 'Parameters okay', 'neutral');
okay=1;

function [DVP, okay] = local_loadDataviewParam(figh, P);
% try to load dataviewparam from file; return void DVP 
okay=1; % optimistic default
if isequal('-', P.Dataviewer), 
    DVP = dataviewparam(); % no viewer active; return void dataviewparam object
elseif isequal('def', lower(P.DataviewerParamfile)),
    DVP = dataviewparam(P.Dataviewer); % default for the current dataviewer 
else, % getting serious: try to load the params from file
    StimType = getGUIdata(figh, 'StimulusType');
    subDir = [StimType '_' P.Dataviewer]; % dir to load & save dataviewparam
    Pfile = strtok(P.DataviewerParamfile,'.');
    DVP = load(dataviewparam(), P.Dataviewer, Pfile, subDir);
    okay = ~isvoid(DVP);
    if ~okay, 
        GUImessage(figh, ...
            {['Dataview parameter file ''' subDir '\' P.DataviewerParamfile ''''], ...
            'not found.'}, 'error', P.handle.DataviewerParamfile);
    end
end    

function  local_edit(Src, nameQ_Pfile, nameQ_Viewer);
% dataviewer param edit button: edit the parameters. Clicking the edit
% button does not depend on the "Check status" of the GUI, so the stimulus
% params are not considered here. Instead, this fcn only uses the 
% Dataviewer and DataviewerParamfile queries and the identity of the
% stimulus GUI.
figh = parentfigh(Src);
GUImessage(figh, ' ');
Q = getGUIdata(figh, 'Query');
StimType = getGUIdata(figh, 'StimulusType');
Viewer = read(Q(nameQ_Viewer));
if isequal('-', Viewer),
    GUImessage(figh, {'No viewer specified;', ...
        'cannot edit its parameters.'}, 'error');
    return;
end
[Pfile, dum, Mess, hFileEdit] = read(Q(nameQ_Pfile), 'coloredit');
if ~isempty(Mess),
    GUImessage(figh, Mess, 'error');
    return;
end
DVP = load(dataviewparam(), Viewer, '?', Pfile); % Pfile is default
if isvoid(DVP),
    GUImessage(figh, 'Editing dataview parameters cancelled.', 'warning');
    return;
end
DVP = edit(DVP);
if isvoid(DVP),
    GUImessage(figh, 'Editing dataview parameters cancelled.', 'warning');
    return;
end
[saved, FN] = save(DVP, '?', Pfile); % prompt user to save ; use fn as default
if saved,
    [dum,fn] = fileparts(FN); % no dir ...
    fn = strtok(fn,'.'); % ... no extension
    GUImessage(figh, ['dataview parameters saved to file ''' fn '''.']);
    set(hFileEdit, 'string', fn);
else,
    GUImessage(figh, 'Saving dataview parameters cancelled.', 'warning');
end




