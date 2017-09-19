function D=get(A);
% address/get - retrieve data from an address
%   D=get(A) retrieves a value previously stored in address A using
%   upload.
%
%   See address/put, address/rm, address/download.

D = A.get();


