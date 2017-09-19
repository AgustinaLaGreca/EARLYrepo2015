function [D, DS, S] = readTKABF(ExpID, RecID, icond, maxNsamStim);
% readTKABF - read ABF file specified by SGSR params
%   [D, DS, s]  = readTKABF(ExpID, RecID, icond);
%    Inputs:
%      ExpID: name of experiment, e.g., 'RG09110'
%      RecID: ID string of SGSR sequence, '3-1-FS'
%      icond: stimulus condition; varied param equals DS.Xval(icond).
%    Outputs:
%       D: ABF data, time aligned and corrected for "missing gain".
%      DS: SGSRdataset
%       s: single entry from TKpool database
%
%    readTKABF(D), where D is a struct, is the same as 
%    readTKABF(D.ExpID, D.RecID, D.icond);
%
%    ReadTKABF adds some extra fields to D compared to readABFdata:
%       ExpID, RecID, icond.
%    readTKABF also applies time alignment (see ABFtimeAlign) and corrects
%    for the "missing gain", i.e., the gain specified directly to the
%    Multiclamp amplifier, which is stored in State files, but not in the
%    header of the raw ABF files.
%    
%      
%
%    See also SGSRdataset, readABFdata, ABFtimeAlign.

if nargin==1 && iscell(ExpID), % unpack cell array
    [ExpID, RecID, icond] = deal(ExpID{:});
    maxNsamStim=inf;
elseif nargin==2 && iscell(ExpID), % 2nd arg is maxNsamStim
    maxNsamStim = RecID;
    [ExpID, RecID, icond] = deal(ExpID{:});
elseif nargin<4, maxNsamStim=inf;
end
%D:\USR\Thomas\data\
persistent TotalPool
if isempty(TotalPool),
    TotalPool = [];
    try,
        TotalPool = [TotalPool, TKpool];
    catch,
        if ~isequal('CLUST', CompuName),
            warning('Problem retrieving TKpool.');
        end
    end
    try,
        TotalPool = [TotalPool, JLpool];
    catch,
        warning('Problem retrieving JLpool.');
    end
end

if isstruct(ExpID), % all args passed in single struct
    [ExpID, RecID, icond] = deal(ExpID.ExpID, ExpID.RecID, ExpID.icond);
end

DS = sgsrdataset(ExpID, RecID);
RecID = DS.SeqID;
S = SGSR_ABF_filter(TotalPool, ExpID, RecID);
% get ABF file info and read file
%BranchDir = strrep(S.Dir, 'D:\USR\Thomas\data\', '');
FFN = fullfile(TKdatadir(S.Experimenter), S.ExpID, S.Pen, S.recID, S.Recs{icond});
MaxNsam = [inf, repmat(maxNsamStim, [1 16])];
D = readABFdata(FFN,1,MaxNsam);
% correct for timing misalignment
D = ABFtimeAlign(D);
% Scale all samples in D using scaling factor 1/S.Rec_Gain. This is 
% typically used for correction of the "missing gain" supplied by the 
% TKpool database (retrieved from state files)
for iad=1:numel(D.AD),
    D.AD(iad).samples = D.AD(iad).samples/S.Rec_Gain;
end



D.ID______________________ = '----------------------------' ;
D.ExpID = ExpID;
D.RecID = RecID;
D.icond = icond;


MyFlag('LastRecordingLoaded', {ExpID, RecID, icond});


