function U = uniqueID(DS);
% Dataset/uniqueID - unique ID of dataset
%    uniqueID(DS) returns the unique ID of dataset DS, which is an integer
%    number  between 1 and 1e15. If DS is an array, uniqueID returns an
%    array having the same size as DS.
%
%    See also dataset.

id = [DS(:).ID];
U = [id.uniqueID];
U = reshape(U, size(DS));









