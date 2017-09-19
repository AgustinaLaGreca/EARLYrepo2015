function [FN, Dir, recID]=ZWOAEfilename(expID, recID, GiveExt, TypeFlag);
% ZWOAEfilename - convert experiment & recording number into ZWOAE datafile name
%
%     ZWOAEfilename(expID, recID) returns the filename of the
%     datafile holding the compacted dataset with number #recID from experiment #expID.
% 
%     ZWOAEfilename(expID, irecID, 0) removes the extension of the 
%     filename. By default, the extension is included.
% 
%     ZWOAEfilename(expID, irecID, [], 'raw') returns file (& path)
%     information for the raw data, rather than the compacted data
% 
%     Examples:
%      ZWOAEfilename(51, 1)   => 'RG07051_zDPs001.ZWcom'
%      ZWOAEfilename(70, 121) => 'RG08070_zDPs121.ZWcom'
%      ZWOAEfilename(70, 121, 1,'raw') => 'RG08070_zDPs121.ZWOAE'
%      ZWOAEfilename(70, 121, 0) => 'RG08070_zDPs121' (3rd arg 0 = no extension)
% 
%     ZWOAEfilename(expID, I) returns a cell array of strings when I is an
%     array.
% 
%     ZWOAEfilename(expID, -14) returns the filename for the system
%     distortion measurement corresponding to recording #14.
% 
%     [FN, Dir]=ZWOAEfilename(expID, recID) also returns the
%     directory Dir where the compacted (or raw) data reside.
% 
%     [FN, Dir, recID]=ZWOAEfilename(expID, inf) returns the filename of
%     the most recently created datafile of the experiment, FN, and its
%     recording index, recID.
% 
%     See also getZWOAEdata, ZWOAEimport.
%

%=== input arg check ===
if nargin<2, recID=[]; end
if nargin<3 || isempty(GiveExt), GiveExt = 1; end % default: provide extention
if nargin<4, TypeFlag=''; end % default: compacted data
if ~isnumeric(expID), 
    error('Input argument ''expID'' must be single integer.');
end

%--- define settings for extensions and subfolders (within zwoaedatadir)---
Comp.Subfolder = ''; %"compacted" data; no extra subfolder
Comp.Ext = '.ZWcom'; %extension for compacted data
Raw.Subfolder = 'RawData'; %subfolder within zwoaedatadir where raw data are saved
Raw.Ext = '.ZWOAE'; %extension for raw data

EXTENSION = Comp.Ext; %default; assume working on compacted data
SUBFOLD = Comp.Subfolder;
if strcmpi(TypeFlag, 'raw'), %no default, working with raw data
    EXTENSION = Raw.Ext;
    SUBFOLD = Raw.Subfolder;
end

%---use recursion for multiple recID---
if numel(recID)>1, 
    for ii=1:numel(recID),
        FN{ii} = ZWOAEfilename(expID, recID(ii), GiveExt, TypeFlag);
    end
    return;
end

% ======single dataset from here=========
SystDist = (recID<1); %when recID<1 --> system distortion
recID = abs(recID);

%decide what year to use in filenames based on experiment number
if expID <= 56,
    YR = '07';
elseif expID <= 85,
    YR = '08';
elseif expID <=188,
    YR = '09';
else,
    YR = '10';
end

expIDStr = dec2base(expID,10,3); % 12 -> '012' etc
FN = ['RG' YR expIDStr];
if expID==0, % TEST
    FN = 'TEST';
end
if SystDist, FN = ['sysdist_for_' FN]; end %system distortion file
Dir = fullfile(ZWOAEdatadir,SUBFOLD,FN);
FN = [FN '_zDPs'];

%----------------------------------------------------------------
%----------------------------------------------------------------
% option to look at frog recordings @ UCLA in 2009 (for Bas Meenderink)
if exist('lookatLA.m','file'),
    if lookAtLA,
        FN = ['LA09' expIDStr];
        if SystDist, FN = ['sysdist_for_' FN]; end %system distortion file
        tmp = ZWOAEdatadir;
        [a b] = fileparts(tmp);
        Dir = fullfile(a,['LA_' b], SUBFOLD, FN);
        FN = [FN '_zDPs'];
    end
end
%----------------------------------------------------------------
%----------------------------------------------------------------

% now incorporate dataset index
if isinf(recID), % % special value of recID: inf==most recent measurement
    [FN, recID]=local_recentst(Dir, EXTENSION);
else, % dataset index is in filename
    FN  = [FN dec2base(recID, 10, 3)];
    FN = fullfilename(FN, '', EXTENSION); % add extension if needed
end

if ~GiveExt, %remove extension when requested
    FN = strtok(FN, '.');
end

%====== LOCAL =======
function [nam, inewest]=local_recentst(Dir, EXT);
[dum Prefix] = fileparts(Dir);
qq = dir(fullfile(Dir,['*' EXT]));
if ~isempty(qq),
    NAME = {qq.name};
    NAM = cellfun(@(s)strtok(s,'.'), NAME, 'UniformOutput',false);
    isep = strfind(NAM, '_zDPs');; % start index - 5 of recID 
    NAM = cellfun(@(s,isep)s(isep+5:end), NAM, isep, 'UniformOutput',false);
    [inewest, ihit] = max(str2double(NAM));
    nam = NAME{ihit};
else,
    inewest = 0;
    nam = '';
end



