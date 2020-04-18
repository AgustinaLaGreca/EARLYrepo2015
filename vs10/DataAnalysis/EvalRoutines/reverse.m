function d = reverse(ds)
%REVERSE returns the dataset such that its spiketrains are reversed in
%        time, i.e. if e.g. the stimulus duration is 1000ms, the old spike
%        times are =t and the new spike times become =1000-t.
% See also: dataset, condition

% Darina Abaffyova 18/04/2020

if isempty(ds) || ~isa(ds, 'dataset') || nargin > 1
    error("First (and only) argument should be a dataset.");
end

% Flip data in the dataset
d = struct(ds);
spt = ds.SpikeTimes;
interval = d.Stimulus.StimParam.interval;

for i = 1:size(spt, 1) % Sequences (rows)
    for j = 1:size(spt, 2) % Repetitions (columns)
        spt{i, j} = abs((spt{i, j}) - interval);
    end
end

d.Data.SpikeTimes = spt;
d = dataset(d, 'convert');