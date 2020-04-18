function d = condition(ds, cond, varargin)
% CONDITION returns the dataset such that it only contains data
%           corresponding to a specified condition
% See also: dataset

% Darina Abaffyova 10/04/2020

if isempty(ds) || ~isa(ds, 'dataset')
    error("First argument should be a dataset.");
end

if isempty(cond) || ischar(cond)
    error("Invalid condition.");
end

if size(cond, 1) > size(cond, 2)
    error("Specify condition as a row vector.");
end

vals = ds.IndepVar.Values;

% If the specified condition doesn't fall in the range of values in the
% dataset
if ~any(cond == vals)
    error("Invalid condition for the dataset.");
end

% Change values in the dataset according to the condition
d = struct(ds);
index = logical(sum(vals == cond, 2));
d.Sizes.NsubRecorded = sum(index);
d.Data.SpikeTimes = ds.SpikeTimes(index, :);
d.Stimulus.IndepVar.Values = ds.IndepVar.Values(index);

% Adapting values in the Stimulus.Special field of the struct
fn = fieldnames(ds.Stimulus.Special);

for i = 1:length(fn)
    f = fn(i);
    temp = eval(['ds.', char(f)]);
    if size(temp) > 1
        d.Stimulus.Special.(char(f)) = temp(index, :); 
    end
end

d = dataset(d, 'convert');