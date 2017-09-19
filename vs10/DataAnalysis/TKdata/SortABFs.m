function [RD, DS, Report] = SortABFs(Experimenter, ExpID);
% SortABFs - compile a list of ABFs and their potential SGSR counterparts
%    SortABFs('TK', '250809') locates all ABF files in the experiment 
%    directory ExpDir = <TKdatadir('TK')>\250809 and lists them in chronological
%    order. Interspersed in the list is the information from the most recent
%    Log files and State files. Log files and State files must be 
%    available in the LOGS and STATEFILES subdirectory of ExpDir.
%
%    [RD, DS, Report] = SortABFs(Experimenter, ExpID) returns the results
%    in variables
%        RD: recursive directory listing (see rdir)
%        DS: struct array containing info per SGGR dataset    
%    Report: char array containing the full report.
%    To save report to file, type 
%       >> textwrite('D:\USR\Jeannette\data\XXXX\ABFreport.log', Report)
%      
%    See also textwrite.

ExpDir = fullfile(TKdatadir(Experimenter), ExpID);
if ~exist(ExpDir,'dir'),
    ['Experiment directory ''' ExpDir ''' not found.']
end

RD = rdir(ExpDir);
RD = sortAccord(RD, [RD.datenum]);


% collect "datasets", i.e. groups of abfs preceeded by a log file
idataset = 0;
StateFile = ''; MultclampState = '???';
LogFile = '????';
Orphans = [];
for ii=1:numel(RD),
    File = RD(ii);
    [dum Name Ext] = fileparts(File.name);
    switch (lower(Ext)),
        case '.state';   % update current state file
            StateFile = File;
            MultclampState = mytextread(File.fullname);
            MultclampState = strvcat(MultclampState{:});
        case '.log', % initialize a new dataset
            %if idataset>2, '--DEBUG--', break; end
            idataset = idataset + 1;
            DS(idataset).StateFile = StateFile;
            DS(idataset).MultclampState = MultclampState;
            DS(idataset).LogFile = File;
            DS(idataset).SGSR_dataset = local_SGSR_DS(File);
            DS(idataset).ABF = []; % no ABFs yet
        case '.abf', % add to current dataset - if any
            if idataset==0,
                Orphans = [Orphans, File];
                continue;
            end
            % store elementary info on ABF file
            abf.name = File.name;
            abf.Nbyte = File.bytes;
            abf.dir = fileparts(File.fullname);
            % read the ABF file & extract basics
            ddd = readABFdata(File.fullname,0); % the 0 means: no samples, just bookkeeping stuff
            abf.Nsweep = ddd.Nsweep;
            abf.sweepDur_ms = ddd.sweepDur_ms;
            abf.mode = ddd.Header.nOperationMode;
            DS(idataset).ABF = [DS(idataset).ABF abf];
        otherwise, % orphans
            Orphans = [Orphans, File];
    end
end

% stats
allABFs = [DS.ABF];
Nabf = numel(allABFs);
NMbyteABF = round(sum([allABFs.Nbyte])/1024^2);
Nds = numel(DS);

% ---report---
Report = ['====SUMMARY of ' ExpID ' ====='];
Report = strvcat(Report, ' ');
Report = strvcat(Report, ['  ' num2str(Nds) ' SGSR data sets.']);
Report = strvcat(Report, ['  ' num2str(Nabf) ' ABF files (' num2str(NMbyteABF) ' Mb).']);
for ii=1:Nds,
    % SGSR params
    Report = strvcat(Report, ' ', '------------------------------------');
    ds = DS(ii);
    Report = strvcat(Report, disp(ds.SGSR_dataset));
    try, Report = strvcat(Report, ['     SPL = ' num2str(ds.SGSR_dataset.SPL) ' dB']); end
    Report = strvcat(Report, ' ');
    % Multiclamp State
    Report = strvcat(Report, ds.MultclampState);
    Report = strvcat(Report, ' ');
    % ABF files
    Na = numel(ds.ABF);
    if Na==0,
        Report = strvcat(Report, '   ####### NO ABF files found ######');
        continue;
    end
    Nrec = ds.SGSR_dataset.Nrec; 
    if Na<Nrec,
        PostFix = [' >>>> ' num2str(Nrec-Na) ' ABFs missing? <<<<<'];
    elseif Na>Nrec,
        PostFix = [' >>>> ' num2str(Na-Nrec) ' ABFs too many? <<<<<'];
    else,
        PostFix = '(regular)';
    end
    Report = strvcat(Report, ['   ' num2str(Na) ' ABF files ' PostFix]);
    Report = strvcat(Report, '   Location(s):');
    ABFdirs = unique({ds.ABF.dir})';
    Report = strvcat(Report, ABFdirs{:});
    for jj=1:numel(ds.ABF),
        abf = ds.ABF(jj);
        LL = ['   ' dec2base(jj,10,2) ':   ' num2str(abf.Nsweep) ' x ' num2str(abf.sweepDur_ms) ' ms;' ];
        LL = [LL ' mode=' num2str(abf.mode) '   ' abf.name];
        Report = strvcat(Report, LL);
    end
end
if ~isempty(Orphans),
    Report = strvcat(Report, '========= ORPHANS & Misc  ==========');
    Ofile = {Orphans.fullname};
    Report = strvcat(Report, Ofile{:});
    Report = strvcat(Report, ' ');
end
Report = strvcat(Report, '========================================');

disp(Report);



%===========================================================
function D = local_SGSR_DS(File);
qq = mytextread(File.fullname);
ww = Words2cell(qq{2}); % second line of log, e.g. 'RG09147 1-1-FS'
try,
    D = sgsrdataset(ww{1}, ww{2});
    D = EmptyDataset(D);
catch,
    %File.fullname
    D.ERROR = ['dataset ' ww{1} ' ' ww{2:end} ' not found;'];
    D.Nrec = 0;
end




