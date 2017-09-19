function DS = SGSRdataset(Exp, recID, Rdir);
% SGSRdataset - constructor for SGSRdataset object. Contains an entire sequence.
%   DS = SGSRdataset('RG09098', '1-1-FS') reads dataset '1-1-FS' of
%   experiment RG09098.
%   You have to export data from SGSR using exportSGSRdata, and copy the
%   resulting directory to rawdata\exportedSGSR. 
%
%   See also SGSRdatadir, SGSRupdata.

if nargin<3,
    Rdir = '';
end
if nargin<1, % void sgsrdataset
    DS = VoidStruct('ID',  'Sizes', 'Data', 'Stimulus', 'Settings');
    DS = class(DS, mfilename);
    return;
end

Ddir = fullfile(SGSRdatadir, ['exported_' , Exp]); % datadir

postFix = strrep(recID, '-', '_'); % postfix of filename
Fname = [Exp '__' postFix '*.SGSRdata'];
FFname = fullfile(Ddir, Fname); % full file name
qq = dir(FFname);
if isempty(qq),
    error(['File ''' FFname ''' not found.']);
elseif numel(qq)>1,
    qq.name
    error(['Multiple files  ''' Fname ''' found.']);
end

FFname = fullfile(Ddir, qq.name);
qq = load(FFname, '-mat');
DS = qq.DS;

DS = class(DS, mfilename);







