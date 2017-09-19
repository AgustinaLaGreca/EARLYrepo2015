function [Z1, Z2]=Z_Pair(rho, Nsam);
%  Z_Pair - two samples of normally distributed variables with given correlation
%    [Z1, Z2]=Z_Pair(rho, Nsam)

if nargin<2, Nsam=1e4; end

Z1 = randn(Nsam,1);
Z2 = randn(Nsam,1);
phi = acos(rho);
Z2 = Z1*cos(phi)+Z2*sin(phi);







