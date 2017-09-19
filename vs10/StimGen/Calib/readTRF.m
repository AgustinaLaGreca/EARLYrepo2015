function TRF = readTRF(FileName, DD, Ext);
% readTRF - read tranfer function from mat file
%   TRF = readTRF(Filename, Dir, Ext) reads a transfer function from file
%   named Filename. Default dir is the 'calibration' subdir of
%   zwoaedatadir. Default extension is '.trf'.
%
%   See also zwoaedatadir, measureTransfer, probeTubeTransfer.


% if TRF is filename, read the file

if nargin<2, DD = ''; end
if nargin<3, Ext = '.trf'; end

if isempty(DD),
    DD = fullfile(ZWOAEdatadir, 'calibration');
end

TRF = fullfilename(FileName, DD, Ext);
if ~exist(TRF,'file'),
    error(['File ''' TRF ''' not found.']);
end
qq=load(TRF);
FN = fieldnames(qq);
if numel(FN)>1,
    error(['File ''' TRF ''' contains multiple variables, not just one transfer function.']);
end
TRF = qq.(FN{1});











