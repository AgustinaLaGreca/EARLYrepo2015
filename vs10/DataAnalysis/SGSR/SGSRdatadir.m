function D = SGSRdatadir;
%  SGSRdatadir - directory where SGSRdata are stored. Hardcoded.

RDD = fileparts(fileparts(EarlyRootDir)); % root dir
if isequal('Cel', CompuName),
    RDD = 'C:\D_Drive';
elseif isequal('CLUST', CompuName),
    RDD = 'D:\';
end
D = fullfile(RDD, 'RawData\exportedSGSR');


