function ds = reverse(ds)
%REVERSE returns the dataset such that its spiketrains are reversed in
%        time, i.e. if e.g. the stimulus duration is 1000ms, the old spike
%        times are =t and the new spike times become =1000-t. Spike times
%        beyond the duration are deleted.
% See also: dataset, dataset/condition, dataset/spiketimes

% Darina Abaffyova 18/04/2020

if isempty(ds) || ~isa(ds, 'dataset') || nargin > 1
    error("First (and only) argument should be a dataset.");
end

PRES = ds.Stim.Presentation;
duration = PRES.PresDur;
spt = spiketimes(ds);

% Use events stored online contained in ds
% for each event, retrieve most recent stimulus onset
ET = digichan(ds, 1);
s = eventtimes(ET);
Onsets = NthFloor(s, PRES.PresOnset(2:end)); % per spike, the most recent onset (ignore 1 = baseline)

% Reverse according to the duration
for i = 1:PRES.Ncond
    dur = duration(i+1); %+1 to account for baseline
    
    for j = 1:PRES.NRep
        spt{i, j} = spt{i, j}(spt{i, j} < dur); % remove spike times beyond stimulus offset
        spt{i, j} = abs(spt{i, j} - dur);
    end
end

SPT = cat(2, spt{:}) + Onsets; % spike times re grand stimulus onset
ds.Data.RX6_digital_1.Data.EventTimes = SPT;

% Remove analog data  if present to avoid confusion (alternative would be
% to reverse these too).
fns = fieldnames(ds.Data);
for f = 1:length(fns)
    if contains(fns(f), 'analog')
        ds.Data = rmfield(ds.Data, fns(f));
    end
end