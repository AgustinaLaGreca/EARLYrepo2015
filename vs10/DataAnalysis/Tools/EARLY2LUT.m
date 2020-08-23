function lut = EARLY2LUT(FN)
%EARLY2LUT get lookup table for EARLY datafile
%   LUT = EARLY2LUT(FN) makes lookup table LUT with all entries from EARLY datafile given by FN
%
%   See also LOG2LUT, EARLYdataset, EARLYdataset

% Darina Abaffyova 20/09/2019

if nargin ~= 1
    error('Wrong number of input arguments.'); 
end

lut = struct('iSeq', {}, 'IDstr', {});

ds = read(dataset, FN, 1); % Is there a better way than having to read the first DS?
list = stimlist(ds.experiment); % This is slowing things down.
for n = 1:numel(list)
    lut(n).iSeq    = n;
    lut(n).IDstr   = list(n).IDstring;
end