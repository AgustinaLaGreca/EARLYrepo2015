function S = OfflineSpikeRecalc(DS);
% Dataset/OfflineSpikeRecalc - redo the calculation of offline spike detection
%    OfflineSpikeRecalc(D) repeats the spike detection & storage for
%    dataset D, using the same function & parameters that was used
%    previously. A warning is given when no spike times were previously
%    stored for D. D may be an array of datasets.
%
%    See also Dataset/OfflineSpikeSave, Dataset/OfflineSpikeRemove, 
%    Dataset/HasOfflineSpikes, Dataset/OfflineSpikeList.

for ii=1:numel(DS),
    ds = DS(ii);
    if HasOfflineSpikes(ds),
        P = OfflineSpikeFullInfo(ds);
        P = P.DetectCall;
        P{2} = ds; % replace dummy by proper dataset
        feval(P{:});
        disp(['Recalculated spike times for ' IDstring(ds,'full')]);
    else,
        warning(['No offline spikes found for ' IDstring(ds,'full') '. Cannot redo the detection.']);
    end
end

















