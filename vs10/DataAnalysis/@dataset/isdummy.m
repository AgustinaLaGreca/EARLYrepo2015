function id=isdummy(D,flag);
% dataset/isdummy - true for dummy dataset
%    isdummy(D) returns true for any dummy dataset.
%
%    isdummy(D,'nostim') returns true only for a dummy dataset that
%    contains no Stim field.
%
%    See also dataset/dummy.

flag = arginDefaults('flag', '');

for ii=1:numel(D),
    id(ii) = ~isvoid(D(ii)) && isfield(D(ii).ID, 'uniqueID') && isequal(0, D(ii).ID.uniqueID);
    if isequal('nostim', flag),
        id(ii) = id(ii) && isequal(struct(),D(ii).Stim);
    end
end
id = reshape(id, size(D));
return;




