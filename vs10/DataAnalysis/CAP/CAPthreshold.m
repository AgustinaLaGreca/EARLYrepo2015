function D = CAPthreshold(expID,iRec)
% CAPthreshold - plot CAP for different recordings of same experiment
%
%   D = CAPthreshold(expID,iRec) plots for experiment expID the CAP data
%   stored in AD1-channel of recordings specified by iRec. The recordings
%   in iRec should be FS and stimulus conditions should only vary in SPL
%   and starting phase (limited to 0 and 0.5). For a given SPL recordings
%   with both starting phases should be included. CAPthreshold first sorts
%   the recordings to their SPL. Second, the mean response for both
%   starting phase 0 and starting phase 0.5 is computed for each SPL. Then,
%   the mean responses of both starting phases are averaged per SPL,
%   removing the CM component and showing the CAP. The user should visually
%   determine the CAP semi-threshold.
%
%   The output struct D has fields:
%       expID       experiment name
%       iRec        recording numbers parsed to the function (iRec)
%       level       sound pressure level for given dataset
%       Nzerophase  number of repetitions used for phase 0
%       Nhalfphase  number of repetitions used for phase 0.5
%       Fcar        stimulus freqeuncies of all conditions
%       data        matrix of CAP data; k-th column is k-th stim. condition
%       dt          dt for time axis
%
%   If for a given SPL multiple recordings with the same starting phase
%   are queried (and one recording has more repetitions than the other), the
%   weighted average is used to get the mean response for that starting
%   phase. The average of the two starting phases is never weighted.
%
%   See also anamean.
%

% Load specified datasets and experiment
cexp = find(experiment(),expID);
DSs = {};
for ii = iRec
    DS = read(dataset,cexp,ii);
    DSs{end+1} = DS;
end

% Check whether datasets are FS
iFS = cellfun(@(DS) strcmp(DS.Stimtype,'FS'),DSs);
if ~all(iFS)
    error('Not all input datasets are of FS type.')
end

% Check whether datasets are similar
simfields = {'StartFreq','StepFreq','StepFreqUnit',...
    'EndFreq','AdjustFreq','FreqTolMode','ModFreq',...
    'ModDepth','ModStartPhase','ModITD','DAC',...
    'BurstDur','RiseDur','FallDur','ITD',...
    'ITDtype','ISI','Baseline',};
nDSs = numel(DSs);
for ii = 1:nDSs
    SPii = structpart(DSs{ii}.Stim,simfields); % struct part index ii
    for jj = ii+1:nDSs % only one-sided similarity matrix necessary, it is more or less mirror symmetric
        SPjj = structpart(DSs{jj}.Stim,simfields); % struct part index jj
        simmat{ii,jj} = structcompare(SPjj,SPii); % cell array representing "similarity matrix"
    end
end
[indI indJ] = find(~cellfun(@isempty,simmat)); % row and column indices that are non-empty and thus need to be returned in error message
if numel(indI) ~=0 % if dissimilar stimulus conditions have been found
    uindJ = unique(indJ); % unique columns from similarity matrix
    mssgstr = {}; % empty error message cell array
    for ii = 1:numel(uindJ)
        logindx = indJ == uindJ(ii); % logical index to current similarity matrix column
        ind = find(logindx); % normal index to current similarity matrix column
        dissimfields = simmat{indI(ind(1)),uindJ(ii)}; % first dissimilar struct for comparison with the rest
        if numel(ind) > 1 % only when more than 1 element of this column is dissimilar
            for jj = 2:numel(ind)
                newfields = simmat{indI(ind(jj)),uindJ(ii)}; % additional dissimilar struct
                dissimfields = structJoin(dissimfields,newfields); % join the two dissimilar structs
            end
        end
        helpvar = cellfun(@(DS) DS.ID.iDataset,{DSs{indI(logindx)}}); % help variable to extract data from cell array
        % Create error message cell array
        mssgstr{end+1} = ['Recording ' ...
            int2str(DSs{uindJ(ii)}.ID.iDataset) ...
            ' is dissimilar from recording(s) ' ...
            strrep(strrep(strrep(int2str(helpvar),'   ','  '),'  ',' '),' ','/') ...
            ' for stimulus condition(s) ' ...
            cell2words(fieldnames(dissimfields).','/')];
    end
    error('Dum:dum',cell2words(mssgstr),'\n') % return multi-line error message
end

% Get SPLs
SPL = cellfun(@(DS) DS.SPL,DSs); % all SPLs
uSPL = unique(SPL); % unique SPLs

% Get starting phase
phi = cellfun(@(DS) DS.WavePhase,DSs); % starting phase for all datasets

% Get number of repetitions
Nrep = cell2mat(cellfun(@(DS) DS.Nrep,DSs,'UniformOutput',false).'); % Nrep for all datasets

% Logical indices to two starting phases
iphi1 = phi == 0;
iphi2 = phi == .5;

D = []; % initialize D
for ii = 1:numel(uSPL)
    cSPL = uSPL(ii); % current SPL
    iSPL = cSPL == SPL; % logical index to current SPLs in recs
    ind1 = all([iphi1;iSPL]); % first combination of SPL & phase
    ind2 = all([iphi2;iSPL]); % second combination of SPL & phase
    % If phase has not been alternated, give a warning and skip the dataset
    if any(ind1)==0 || any(ind2)==0
        warning(['No 0-cycle & 0.5-cycle start phase available for intensity of ' num2str(cSPL) ' dB SPL (iRec ' int2str(iRec(any([ind1;ind2]))) '). SPL skipped.']) %#ok<WNTAG>
        continue
    end
    cDSs1 = {DSs{ind1}}; % current dataset(s) 1
    cDSs2 = {DSs{ind2}}; % current dataset(s) 2
    temp1 = Nrep(ind1,:); % number of repetitions per condition for dataset(s) 1
    for jj = 1:numel(cDSs1)
        [D1(:,:,jj) dt] = anamean(cDSs1{jj},1); % call anamean for dataset(s) 1 and also get dt
        weights1(:,:,jj) = SameSize(temp1(jj,:),D1(:,:,jj)); % get weights for dataset(s) 1
    end
    temp2 = Nrep(ind2,:); % number of repetitions per condition for dataset(s) 2
    for jj = 1:numel(cDSs2)
        D2(:,:,jj) = anamean(cDSs2{jj},1); % call anamean for dataset(s) 2
        weights2(:,:,jj) = SameSize(temp2(jj,:),D2(:,:,jj)); % get weights for dataset(s) 2
    end
    % Compute mean using weights obtained by number of repetitions
    D1sum = sum(D1.*weights1,3);
    D2sum = sum(D2.*weights2,3);
    sumweights1 = sum(weights1,3);
    sumweights2 = sum(weights2,3);
    Dmean(:,:,1) = D1sum./sumweights1;
    Dmean(:,:,2) = D2sum./sumweights2;
    % Put resulting data in struct
    D(end+1).expID = expID;
    D(end).iRec = iRec(any([ind1;ind2]));
    D(end).level = cSPL;
    D(end).Nzerophase = sumweights1(1,:);
    D(end).Nhalfphase = sumweights2(1,:);
    D(end).Fcar = cDSs1{1}.Fcar;
    D(end).data = mean(Dmean,3);
    D(end).dt = dt;
    % Clear some variables
    clear Dmean weights1 weights2 D1 D2
end

% If no level with 0 and 0.5 start phase has been found, return
if numel(D) == 0, return, end

% Call plot command
local_plot(D)


function local_plot(D)
% LOCAL_PLOT - plots CAP-data from output struct
%

% Figure "dimensions"
nrows = numel(D); % number of rows
ncols = 4; % number of columns
nCond = numel(D(1).Fcar); % number of conditions
nfigs = ceil(nCond/4); % number of figures
nsubp = nrows*ncols; % number of subplots

% Plot CAP-data struct condition by condition
for ii = 1:nCond
    cfig = floor((ii-1)/ncols)+1; % current figure
    ccol = rem((ii-1),ncols)+1; % current column
    figure(cfig) % activate figure
    for jj = 1:nrows
        subplot(nrows,ncols,ccol+(jj-1)*ncols) % select correct subplot
        dplot(D(jj).dt,D(jj).data(:,ii),[ploco(jj) '-']) % plot
        legend({[num2str(D(jj).level) ' dB SPL; ' int2str(D(jj).Nzerophase(ii)) '/' int2str(D(jj).Nhalfphase(ii))]},'Location','SouthEast') % add legend
        if jj == 1
            title([num2str(D(jj).Fcar(ii)/1e3,2) ' kHz']) % add title when necessary
        end
    end
end

% Set labels
for ii = 1:nfigs
    figure(ii)
    axes
    set(gca,'Color','none','XColor',[.8 .8 .8],'YColor',[.8 .8 .8])
    xlabel('Time (ms)','Color','k')
    ylabel('Response (V)','Color','k')
    ch = get(gcf,'Children');
    uistack(ch(1),'bottom')
end



















