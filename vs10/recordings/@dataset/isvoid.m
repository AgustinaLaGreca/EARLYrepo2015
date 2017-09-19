function v=isvoid(D);
% Dataset/isvoid - true for void dataset.
%   isvoid(D) returns TRUE if D is a void dataset, FALSE otherwise. A void
%   dataset is a dataset returned by the dataset constructor called without
%   any input arguments. If D is an array of datasets, isvoid returns a
%   corresponding logical array.
%
%   See Dataset.

for ii=1:numel(D),
    v(ii) = isempty(D(ii).ID);
end
v = reshape(v, size(D));
