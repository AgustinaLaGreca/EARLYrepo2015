function D=figh(A);
% address/figh - retrieve graphics handle from an address
%   D=figh(A) retrieves the graphics handle of teh graphics object where A 
%   is stored. If A is not stored in a graphics object, figh returns [].
%
%   See address/put, address/rm, address/download.

D = A.figh;


