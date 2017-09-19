function S = OfflineSpikeInfo(DS);
% Dataset/OfflineSpikeInfo - get additional info saved with offline spikes 
%    S = OfflineSpikeInfo(D) returns the struct S previously saved by 
%    OfflineSpikeSave. An empty struct is returned when no offline spikes 
%    times were saved for dataset D.
%
%    See also Dataset/OfflineSpikeSave, Dataset/OfflineSpikeRemove, 
%    Dataset/HasOfflineSpikes, Dataset/OfflineSpikeList,


if ~HasOfflineSpikes(DS),
    error(['No offline spiketimes found for dataset ' IDstring(DS,'full') '.']);
end
S = OfflineSpikeFullInfo(DS);
S = S.AddInfo;

















