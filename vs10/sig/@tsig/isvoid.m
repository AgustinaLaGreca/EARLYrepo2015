function IV = isvoid(T)
% tsig/isvoid - true for void tsig object
%   isvoid(T) tests if T is a void tsig object, i.e., he result of a call
%   to tsig with no input arguments.
%
%   See tsig.

IV = isempty(T.Fsam); 





