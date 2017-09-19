function [Sall, Date] = TKpool(flag, Experimenter);
%  TKpool - database of TK's ABF data
%    S = TKpool returns struct array S containing info that uniquely links
%    Thomas' ABF files with their SGSR counterparts.
%    S = TKpool('return') is the same thing.
%
%    TKpool('update') attempts to update the pooled data from KULAK.
%
%    TKpool('list') lists the experiments contained in the pool.
%
%    See also JLpool, ABFupdata, SGSRupdata, compile_all_SGSR_ABFs, zipwork.

if nargin<1, flag='return'; end
if nargin<2, Experimenter = 'TK'; end

persistent SSS DDD EEE

more off
[flag, Mess] = keywordMatch(flag, {'list' 'update' 'return'}, 'flag');
error(Mess);
switch flag,
    case 'list', % list elementary info
        [qq, dd] = TKpool('return', Experimenter);
        [SGSRname, iun] = unique({qq.SGSRExpID});
        ABFname = {qq(iun).ExpID};
        sep = repmat({' = '}, 1,numel(iun));
        AllNames = [SGSRname; sep; ABFname]';
        disp(['--- TKpool was last updated ' dd]);
        disp(cell2chararray(AllNames));
        disp(['---']);
        local_listNew(ABFname);
    case 'update', % delegate to compile_all_SGSR_ABFs
        local_update(Experimenter);
        [SSS, DDD, EEE] = deal([]); % erase any previous persistent value
        clear readTKABF; % clear persistent pool in readTKABF
    case 'return', % return the pool 
        if ~isempty(SSS) && isequal(EEE, Experimenter), % return persistent  versions
            Sall = SSS;
            Date = DDD;
        else, % load from file
            PoolFile = fullfile(TKdatadir(Experimenter), 'SGSR_ABF_link', 'All_Exps.SGSR_ABF');
            Sall = load(PoolFile, '-mat'); Sall = Sall.S;
            Date = dir(PoolFile); Date = Date.date;
            SSS = Sall;
            DDD = Date;
            EEE = Experimenter;
        end
        Sall = localRemoveFails(Sall); % remove datasets whose bookkeeping is corrupted due to communication errors
end

%===================

function CommonPool = local_commonpool,
if isequal('KULAK', CompuName), % special case: no network drive, but local dir
    CommonPool = 'D:\MatlabProgs\TKpool_info\all.SGSR_ABFs';
else, % from kulak over network
    CommonPool = 'K:\TKpool_info\all.SGSR_ABFs';
end

function local_listNew(ABFname);
for ii=1:numel(ABFname),
    an = ABFname{ii};
    fan = fullfile(TKdatadir,an);
    if ~exist(fan, 'dir'),
        disp(['-------->>> ', an, ' not uploaded']);
    end
end


function S = local_update(Experimenter);
% local_update - combine all SGSR_ABF files and save.
%   compile_all_SGSR_ABFs(Experimenter) looks for all experiments in
%   TKdatadir(Experimenter) and compiles a database linking the ABF files
%   to their SGSR counterparts by calling SGSR_ABF_link.
% 
%  See also SGSR_ABF_link, ABFupdata, SGSRupdata.
if nargin<1, Experimenter='TK'; end % Thomas the pioneer.
Ddir = TKdatadir(Experimenter);
% find out which experiments are in Ddir.
qq = dir(Ddir);
qq = qq([qq.isdir]);
qq = qq(3:end); % remove .\ and ..\
S = [];
for ii=1:numel(qq),
    xp = qq(ii).name
    if exist(fullfile(Ddir, xp, 'DontPoolMe.txt')),
        disp(['*** Skipped ' xp ' because of DontPoolMe.txt tag file. ***']);
    elseif ~isequal('SGSR_ABF_LINK', upper(xp)) && ~isequal('STIMULUSLISTS', upper(xp)),
        [s, Mess] = SGSR_ABF_link(xp);
        warning(Mess);
        if ~isempty(s),
            if isempty(S), S=s; else, S=[S,s]; end
        end
    end
end
% combine all .SGSR_ABF files in pool
PoolFile = fullfile(TKdatadir(Experimenter), 'SGSR_ABF_link' , 'All_Exps.SGSR_ABF');
save(PoolFile, 'S');

function  S = localRemoveFails(S); 
% remove datasets whose bookkeeping is corrupted due to communication
% errors of the hybrid setup (missed synch pulses, SGSR crashes, etc)
iselect = strmatch('RG10197', {S.ExpID});
ibug = strmatch('20:28:00', {S(iselect).Timestamp_first_abf});
S(iselect(ibug)) = [];
















