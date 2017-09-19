function p = get(P, pnam);
% dataviewparam/get - get parameters from dataviewparam objects.
%   get(P) returns the Param struct in P, holding the parameter values.
%   get(P, 'Foo') returns parameter Foo.

p = P.Param;
if nargin>1,
    p = p.(pnam);
end




