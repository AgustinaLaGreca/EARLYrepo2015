function X=numel(H);
% hoard/numel - number of elements of data of hoard object
%    numel(H) returns the number of elements of the data of hoard object H.
%
%    NOTE: this includes both samples stored in memory and on disk!
%
%    See also hoard/set, hoard/cat.

UD = get(H.h, 'userdata');
X=numel(UD.data) + sum(UD.NsamSaved);







