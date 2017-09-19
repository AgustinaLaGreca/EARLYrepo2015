function f=tempfileprefix(H,X, Dim);
% hoard/tempfileprefix - file name prefix for temporary hoard data storage
%    tempfileprefix(H) returns full path and filename prefix of temporary
%    hoard files. These are used by H to dump long arrays.
%
%    See also hoard/cat, hoard/cleanup.

UD = get(H.h, 'userdata');
f = [UD.Filename '_'];



