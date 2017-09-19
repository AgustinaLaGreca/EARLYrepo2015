function INL = isNumOrLogical(X);
% isNumOrLogical - true for numeric or logical arrays
%    isNumOrLogical(X) returns logical true (1) if X is a numeric or a
%    logical value, and logical false (0) otherwise.
%
%    See also ISNUMERIC, ISLOGICAL.

INL = isnumeric(X) || islogical(X);

