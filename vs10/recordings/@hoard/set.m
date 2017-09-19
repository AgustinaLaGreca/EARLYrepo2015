function set(H,X);
% hoard/set - set data of hoard object
%    set(H,X) sets data of hoard H to X, overwriting any previous value.
%
%    Note that calling set(H,{}) makes H's data a cell array, affecting 
%    subsequent cat(H,..) calls.
%
%    See also hoard/set, hoard/cat.

if nargin<3, Chunk = []; end
UD = get(H.h, 'userdata');
UD.data = X;
set(H.h, 'userdata', UD);
N = numel(UD);







