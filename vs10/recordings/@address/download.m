function D=download(A);
% address/download - retrieve data from an address
%   D=download(A) retrieves a value previously stored in address A using
%   upload. It is completely equivalent to D=get(A).
%
%   See address/put, address/rm.

D = A.get();


