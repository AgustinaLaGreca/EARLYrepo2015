function FN = OfflineSpikeFilename(DS, CreateDir);
% Dataset/OfflineSpikeFilename - filename for offilne spike times
%    FN = OfflineSpikeFilename(DS, CreateDir)
%    returns full filename for offline spike times of dataset DS. 
%    If CreateDir==true, the experiment subdir for this file is created if
%    needed.
%
%    See also Dataset/OfflineSpikeRemove, Dataset/HasOfflineSpikes,
%    Dataset/OfflineSpikeInfo, Dataset/OfflineSpikeList.

[CreateDir] = arginDefaults('CreateDir', false);


% Each exp has its own subdir. Provide it if needed.
SubDir = 'Early_spiketimes';
Sdir = fullfile(processed_datadir, SubDir, expname(DS));
if CreateDir && ~exist(Sdir,'dir'),
    [ok, mess] = mkdir(fullfile(processed_datadir, SubDir), expname(DS));
    error(mess);
end

FN = fullfile(Sdir, ['DS', zeropadint2str(irec(DS),4), '.offspt']);


















