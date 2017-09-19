function [S, Mess] = SGSR_ABF_link(ExpID)
% SGSR_ABF_link - link SGSR & ABF data using hybrid_bookkeeping
%   S=SGSR_ABF_link('160409') returns struct array S summarizing all data in 
%       the experiment folder <TKdatadir>\160409.
%   S=SGSR_ABF_link('RG09163') returns struct array S summarizing all data in 
%       the experiment folder <JLdatadir>\RG09163.
%   Struct array S is also saved in file XXX.SGSR_ABF in the SGSR_ABF_link
%   subfolder of the experiment folder.
%   
%   Example:
%   
%   SGSR_ABF_link('RG09163')
%   ans = 
%   1x29 struct array with fields:
%       ExpID
%       Experimenter
%       Pen
%       recID
%       irec
%       Dir
%       Recs
%       qq______________________
%       Date_of_recording
%       Timestamp_first_abf
%       Timestamp_last_abf
%       d_____
%       SGSRExpID
%       SGSRidentifier
%       Stim_Property_Unit
%       Stim_Property_Values
%       Rec_Quality
%       Rec_Base_unit
%       Rec_Gain
%
%   See also TKdatadir, hybrid_bookkeeping.

Mess = '';
DdirTK = TKdatadir('TK');
DdirJL = TKdatadir('JL');
% find out who's calling
if exist(fullfile(DdirTK, ExpID), 'dir'),
    Ddir = DdirTK;
    UnitPrefix = 'pen*';
    RecordPrefix = 'record*';
    Experimenter = 'TK';
elseif exist(fullfile(DdirJL, ExpID), 'dir'),
    Ddir = DdirJL;
    UnitPrefix = 'unit*';
    RecordPrefix = 'record*';
    Experimenter = 'JL';
else, 
    error(['Cannot find ABF data of experiment ''' ExpID '''.']);
end

RootDir = fullfile(Ddir, ExpID);
Plist = dir(fullfile(RootDir, UnitPrefix));
S = [];
for ipen=1:numel(Plist),
    PenName = Plist(ipen).name;
    Rlist  = dir(fullfile(RootDir, PenName, RecordPrefix));
    for irec=1:numel(Rlist),
        RecName = Rlist(irec).name;
        FullName = fullfile(RootDir, PenName, RecName);
        % get indiv recordings & sort according to date
        qq = dir(fullfile(FullName,'*.abf'));
        [dum, isort]  = sort([qq.datenum]); 
        qq = qq(isort);
        % store info on ABF recording
        s1.ExpID = ExpID;
        s1.Experimenter = Experimenter;
        s1.Pen = PenName;
        s1.recID = RecName;
        s1.irec = local_TrailNum(RecName);
        s1.Dir = FullName;
        s1.Recs = {qq.name};
        s1.qq______________________ = '__________________';
        % combine with SGSR info 
        s2 = hybrid_bookkeeping(FullName);
        if ~isempty(s2),
            s = structJoin(s1,s2);
            S = [S, s];
        end
    end
end

if isempty(S),
    Mess = ['SGSR-ABF link failed for ''' ExpID '''.'];
    return;
end

SGSRid = S(1).SGSRExpID; % SGSR name of experiment
% save the info 
%SaveDir = fileparts(which(mfilename));
SaveDir = fullfile(Ddir,'SGSR_ABF_link');
if ~exist(SaveDir,'dir'), % create it
    [dum, Mess] = mkdir(Ddir,'SGSR_ABF_link');
    error(Mess);
end
SaveName = fullfile(SaveDir, [SGSRid '.SGSR_ABF']);
save(SaveName, 'S', '-mat');

%====================
function x=local_TrailNum(str);
qisNum = betwixt(str, '0'-1, '9'+1);
ilast = find(~qisNum, 1, 'last'); % index of last non-numeric char in str
x = str2num(str(ilast+1:end));











