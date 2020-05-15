function ds = reverse(ds)
%REVERSE returns the dataset such that its spiketrains are reversed in
%        time, i.e. if e.g. the stimulus duration is 1000ms, the old spike
%        times are =t and the new spike times become =1000-t. Spike times
%        beyond the duration are left as-is.
% See also: dataset, dataset/condition, dataset/spiketimes

% Darina Abaffyova 18/04/2020

if isempty(ds) || ~isa(ds, 'dataset') || nargin > 1
    error("First (and only) argument should be a dataset.");
end

PRES = ds.Stim.Presentation;
BurstDur = ds.Stim.BurstDur;

% Use events stored online contained in ds
% for each event, retrieve most recent stimulus onset
ET = digichan(ds, 1);
spt = eventtimes(ET);
[Onsets, ipres] = NthFloor(spt, PRES.PresOnset(2:end)); % per spike, the most recent onset (ignore 1 = baseline)
spt = spt - Onsets; % spike times re "their own" stim onsets

% Reverse according to the duration
for i = 1:PRES.Ncond*PRES.NRep
    spt(ipres == i & spt < BurstDur) = abs(spt(ipres == i & spt < BurstDur) - BurstDur);
end
spt = spt + Onsets; % spike times re grand stimulus onset
spt(spt > BurstDur + Onsets) = [];
ds.Data.RX6_digital_1.Data.EventTimes = spt;

% Remove analog data  if present to avoid confusion (alternative would be
% to reverse these too).
fns = fieldnames(ds.Data);
for f = 1:length(fns)
    if contains(fns(f), 'analog')
        ds.Data = rmfield(ds.Data, fns(f));
    end
end