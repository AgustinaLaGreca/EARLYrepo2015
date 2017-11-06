function headerStr = getDefaultHeaderText(ds)
% GETDEFAULTHEADERTEXT returns the default text for a HeaderObject

%% first column
filename = upper(ds.ID.Experiment.ID.Name);
fileformat = upper('EarlyDS');
seqid = sprintf('%s (%d)', [num2str(ds.ID.iCell) '-' num2str(ds.ID.iRecOfCell) '-' ds.StimType], ds.ID.iDataset);
stimtype = upper(ds.StimType);
% rectime = datestr(datenum(ds.ID.created(1:2), ds.ID.created(4:6), ds.ID.created(8:11)), 1);
month = [];
switch ds.ID.created(4:6)
    case 'Jan'
        month = 1;
    case 'Feb'
        month =2;
    case 'Mar'
        month= 3;
    case 'Apr'
        month = 4;
    case 'May'
        month = 5;
    case 'Jun'
        month = 6;
    case 'Jul'
        month = 7;
    case 'Aug'
        month = 8
    case 'Sep'
        month = 9;
    case 'Oct'
        month = 10
    case 'Nov'
        month = 11;
    case 'Dec'
        month = 12;
end

rectime = datestr(datenum(str2num(ds.ID.created(8:11)), month, str2num(ds.ID.created(1:2))), 1);

firstCol = sprintf([' FileName \n FileFormat \n SeqID (iSeq) \n StimType ' ...
    '\n RecTime | %s \n %s \n %s \n %s \n %s '], ...
    filename, fileformat, seqid, stimtype, rectime);

%% second column

    labels = {' Channels \n', ' NSubSeqs/NReps \n', ' IndepVar \n', ...
        ' IndepRange \n', ' SPL | %s \n %s \n %s \n %s \n %s '};
    values = {sprintf('%d/%d x %d', ds.Stim.Presentation.Ncond, ds.Stim.Presentation.Ncond,ds.Stim.Presentation.Nrep), ...
        ds.Stim.Presentation.X.ParName, indepVar2Str(ds.Stim.Presentation.X), spl2Str(ds)};

secondCol = sprintf([labels{:}], values{:});

%% third column
%To be sure frequency conventions are equal between EDF and
%SGSR/Pharmington datasets ...
S = GetFreq(ds);
%Displaying only the appropriate beat frequency ...

nrec = ds.Stim.Presentation.Ncond;
if all(isnan(S.BeatFreq(1:nrec))) && all(isnan(S.BeatModFreq(1:nrec)))
    BeatFreq = NaN;
elseif ~all(isnan(S.BeatModFreq(1:nrec)))
    %Beat on modulation frequency has priority over carrier frequency ...
    BeatFreq = S.BeatModFreq;
else
    BeatFreq = S.BeatFreq;
end
thirdCol = sprintf([' BurstDur \n RepDur \n CarFreq \n ModFreq \n BeatFreq |' ...
    ' %s \n %s \n %s \n %s \n %s '], ...
    num2str(ds.Stim.BurstDur),num2str(ds.Stim.ISI), ...
    freqVar2Str(S.CarFreq(1:nrec, :)), freqVar2Str(S.ModFreq(1:nrec, :)), ...
    freqVar2Str(BeatFreq));

%% putting it all together
headerStr = [firstCol '|' secondCol '|' thirdCol];
