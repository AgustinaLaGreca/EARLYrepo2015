function S = OfflineSpikeFullInfo(DS);
% Dataset/OfflineSpikeFullInfo - retrieve full info of offline spike detection
%    S = OfflineSpikeInfo(D) returns the full info saved by 
%    OfflineSpikeSave. An empty struct is returned when no offline spikes 
%    times were saved for dataset D.
%
%    See also Dataset/OfflineSpikeSave, Dataset/OfflineSpikeRemove, 
%    Dataset/HasOfflineSpikes, Dataset/OfflineSpikeList,
%    Dataset/OfflineSpikeInfo.

P = struct([]);
FN = OfflineSpikeFilename(DS); 
if exist(FN, 'file'),
    S = load(FN, '-mat');
end


















