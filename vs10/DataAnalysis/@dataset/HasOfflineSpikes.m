function H = HasOfflineSpikes(DS);
% Dataset/HasOfflineSpikes - true if offline spikes were saved for dataset
%    HasOfflineSpikes(D) returns True if offline spike times were saved for
%    dataset D. This hinges on the existence of the file in which these
%    spike times are saved.
%
%    See also Dataset/OfflineSpikeSave, Dataset/OfflineSpikeRemove,
%    Dataset/OfflineSpikeFilename, Dataset/OfflineSpikeList.

H = logical(exist(OfflineSpikeFilename(DS), 'file'));



















