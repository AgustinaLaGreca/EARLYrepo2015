function D = read(DS, ExpName, iRec); %#ok<INUSL,STOUT>
% Dataset/read - read Dataset from disk.
%    DS = read(dataset, XP, iRec) where XP is an Experiment object, reads
%    dataset # iRec from experiment XP. Experiments may also be specified
%    by their name:
%    DS = read(dataset, 'RG09211', 18) reads dataset # 18 from experiment RG09211.
%
%    DS = read(dataset, XP, inf) reads most recently saved dataset.
%
%    DS = read(dataset, XP, I), where I is an array of indices, reads
%    multiple datasets into dataset array DS.
%
%    DS = read(dataset, XP, 0), reads all XP datasets into dataset array DS.
%
%    DS = read(DS, 0) "updates" DS by appending all datasets of
%    experiment(DS) that are not yet in DS.
%
%    DS = read(DS, N) appends only datasets upto #N.
%
%
%    DS = read(dataset, ID), where ID is the struct returned by a call to
%    dataset/ID, returns the original dataset from which ID was extracted.
%
%    See Dataset.

% handle all multiple reads first by recursion
if ~isvoid(DS), % "update", i.e., append to existing DS  (see help)
    % error('reading dataset arrays not yet implemented');
    if nargin>2, error('Too many input args for non-void dataset input. See help.'); end
    D = DS;
    iRec = ExpName;
    ExpName = name(DS(1).ID.Experiment);
    EE = find(experiment, ExpName); % up-to-date version of experiment
    Nrec = status(EE,'Ndataset'); % # datasets sofar recorded for exp
    if isequal(0, iRec), % add all missing
        iRec = setdiff(1:Nrec, irec(DS)); % indices of recorded datasets not in DS yet
    elseif isscalar(iRec), % N means 1:N (see help)
        iRec = setdiff(1:iRec, irec(DS)); %
    end
    for ii=1:numel(iRec),
        D(end+1) = read(dataset, ExpName, iRec(ii));
    end
    return
elseif numel(iRec)>1, % fully specified index array 
    for ii=1:numel(iRec),
        D(ii) = read(dataset, ExpName, iRec(ii));
    end
    return;
elseif isequal(0, iRec), % by convention, zero means "all" datasets
    EE = find(experiment, ExpName); % up-to-date version of experiment
    Nrec = status(EE,'Ndataset');
    for ii=1:Nrec,
        D(ii) = read(dataset, ExpName, ii);
    end
    return
end


%===========single dataset from here =================
if isa(ExpName, 'experiment'),
    if isvoid(ExpName),
        error('Cannot read dataset from a void experiment.');
    end
    ExpName = name(ExpName);
elseif isstruct(ExpName),
    [ExpName, iRec] = deal(name(ExpName.Experiment), ExpName.iDataset);
end

Ddir = locate(experiment, ExpName);
if isinf(iRec), % last one
    EE = find(experiment, ExpName);
    iRec = status(EE,'Ndataset');
end
FileName = [ExpName '_DS' zeropadint2str(iRec,5) '.EarlyDS'];
FFN = fullfile(Ddir, FileName);
%if local_isOldgerbil(ExpName), % suppress warning on experiment object def, which is fixed below
    warning('off', 'MATLAB:elementsNowStruc');
%end
load(FFN, 'D', '-mat'); 
warning('on', 'MATLAB:elementsNowStruc');

D = local_fix(D);

%=========================================
function D = local_fix(D);
% fix known bookkeeping bugs
if isequal('MASK', D.Stim.StimType) && ~isfield(D.Stim, 'DurOkay'), % fix bug in burstdur storage
    D.Stim.Duration = D.Stim.Duration.*4.^(-D.Stim.Warp); 
end
if isequal('MASK', D.Stim.StimType) && size(D.Stim.Duration,2)>1, % fix bug in "fixed" burstdur storage
    WarpFactor = 2.^D.Stim.Warp;
    D.Stim.Duration = max([D.Stim.MnoiseBurstDur (D.Stim.ToneOnset+D.Stim.ToneDur)])./WarpFactor;
end
% added new ReducedStorage field to experiment object; add it to ensure
% backward compatibility
D.Stim.Presentation = stimpresentx(struct(D).Stim.Presentation);
% added new Sys3setup field to experiment object; add it to ensure backward compatibility
D.ID.Experiment = experiment(struct(D).ID.Experiment);
D.Stim.Experiment = experiment(struct(D).Stim.Experiment);

function IOG = local_isOldgerbil(ExpName)
IOG = false;
if numel(ExpName)>2 && isequal('RG', upper(ExpName(1:2))) && str2num(ExpName(3:end))<=11306,
    IOG = true;
end

