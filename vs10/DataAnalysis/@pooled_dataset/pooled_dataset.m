function P = pooled_dataset(DS);
% pooled_dataset - construct pooled_dataset object
%   P = pooled_dataset(DS) constructs a pool of datasets from a dataset 
%   array DS.

if nargin<1, DS = dataset(); end
P = struct('DS', DS);
P = class(P, mfilename, dataset('pooled_dataset'));


