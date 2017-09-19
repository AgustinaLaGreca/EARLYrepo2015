function okay=save(D,varargin);
% Dataset/save - save Dataset object.
%
%   See Dataset.


eval(IamAt);

if ~isequal('complete', D.Status), D = download(D); end;

if isfield(D.ID,'donotsave')
   D=rmfield(D.ID,'donotsave');
   return;
end

% interrupted data may or may not be saved - ask user
if isequal('interrupted', D.Status),
    YN = questdlg('Save interrupted recording?', 'User interrupt during recording.', 'Yes', 'No', 'Yes');
    if isequal('No', YN), % clean up any temp data on disk & beat it
        cleanup(D);
        GUImessage(gcg, {'Recording cancelled;',  'no data saved.'});
        okay=true; 
        return; 
    end
elseif ~isequal('complete', D.Status),
    error(['Cannot save Dataset object whose status is ''' D.Status '''.']);
end

% get status of current exp 
if isempty(varargin)
    EE = D.ID.Experiment;
    StE = status(EE); % status of experiment prior to this save;
    iDataset = StE.Ndataset+1;
    D.ID.priorExpStatus = StE;
    D.ID.iDataset = iDataset;
    D.ID.iCell = max(1, StE.iCell);
    D.ID.iRecOfCell = StE.iRecOfCell+1;
    D.ID.iPen = StE.iPen;
    D.ID.PenDepth = StE.PenDepth;
    D.ID.contributors = [D.ID.contributors, fullfile(mfilename, mfilename('class'))]; % sign your work
elseif strcmp(varargin{1},'addcomment')
    EE = D.ID.Experiment;
    StE = status(EE); % status of experiment prior to this save;
    iDataset = StE.Ndataset;
    D.ID.priorExpStatus = StE;
    D.ID.iDataset = iDataset;
    D.ID.iCell = max(1, StE.iCell);
    D.ID.iRecOfCell = StE.iRecOfCell;
    D.ID.iPen = StE.iPen;
    D.ID.PenDepth = StE.PenDepth;
    D.ID.contributors = [D.ID.contributors, fullfile(mfilename, mfilename('class'))]; % sign your work
end
% some types of data, e.g. ADC data, may have been (partially) dumped to 
% disk. This requires post-handling: storing the data under appropriate 
% names, removing remaining samples from D, and providing D with the 
% information on the disk storage. This is delegated to Save methods of the
% data themselves.
DFNS = fieldnames(D.Data);
for ii=1:numel(DFNS),
    dfn = DFNS{ii};
    if (~strcmpi(dfn,'Thr')&...
        ~strcmpi(dfn,'Freq_masker')&...
        ~strcmpi(dfn,'dB_masker')&...
        ~strcmpi(dfn,'thrs')&...
        ~strcmpi(dfn,'thr_probe')&...
        ~strcmpi(dfn,'cf')&...
        ~strcmpi(dfn,'thr_data')&...
        ~strcmpi(dfn,'SR'))
        D.Data.(dfn) = save(D.Data.(dfn), D, dfn);
    end
end
% save D in MAT format
temp_experimenter = experimenter(D.ID.Experiment);
D.ID.Experiment= setExperimenter(D.ID.Experiment,'Exp');
warning('off','MATLAB:unassignedOutputs');
save([fileprefix(D) '.EarlyDS'], 'D', '-mat');
warning('on','MATLAB:unassignedOutputs');
D.ID.Experiment = setExperimenter(D.ID.Experiment,temp_experimenter);
% update log file and experiment status
id = D.ID;
IDstr = ['Dataset # ' num2str(id.iDataset) ':  ' , ...
    IDstring(D), '     ', id.modified];
addtolog(EE, IDstr);
AllSaved = status(EE, 'AllSaved'); % AllSaved(k,:) = [icell, iofcell] of k-th set saved
EE = state(EE, 'Recording');
StimType = D.Stim.StimType; iCell = D.ID.iCell; iRecOfCell = D.ID.iRecOfCell;
status(EE, 'Ndataset', id.iDataset, ...
    'IDlastSaved', IDstring(D), ...
    'AllSaved', [AllSaved; CollectInStruct(iDataset, iCell,  iRecOfCell, StimType)], ...
    'iRecOfCell', id.iRecOfCell, ...
    'iCell', id.iCell);

GUImessage(gcg, ['Dataset '  num2str(id.iDataset)  ' saved [' IDstring(D) '].']);
okay = true;


