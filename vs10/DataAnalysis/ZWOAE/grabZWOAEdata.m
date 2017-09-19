function N = grabZWOAEdata(iGerbil);
% grabZWOAEdata - download raw data of ZWOAE experiment
%   grabZWOAEdata(iGerbil) copies compacted data from YXO to local rawdata
%   directory. The YXO datadir is supposed to be a local mapped 
%   network drive M:\.

% get file list of avaialable data
FileFilter = ZWOAEfilename(iGerbil,9999999);
FileFilter = strrep(FileFilter, '9999999', '*');
subdir = strtok(FileFilter, '_');
SourceDir = fullfile('M:\ZWOAE', subdir);
% SourceDir = fullfile('M:\ZWOAE', subdir, 'compact');
Sdata = dir(SourceDir);
Sdata = Sdata(~[Sdata.isdir]);
if isempty(Sdata), % escape here - don't create local rawdata dir for nothing
    warning(['No ZWOAE data for gerbil ' subdir '.']);
    N = 0;
    return;
end
[dum, isort] = sort([Sdata.datenum]); % oldest first
Sdata = Sdata(isort);

% find out what we have
DestDir = fullfile(ZWOAEdatadir, subdir);
if ~exist(DestDir, 'dir'), % create it
    [okay, Mess] = mkdir(ZWOAEdatadir, subdir);
    error(Mess);
end

% copy only data that we don't have yet
for ii=1:numel(Sdata),
    Src = fullfile(SourceDir, Sdata(ii).name);
    Dest = fullfile(DestDir, Sdata(ii).name);
    if ~exist(Dest,'file'),
        copyfile(Src, DestDir);
    end
end






