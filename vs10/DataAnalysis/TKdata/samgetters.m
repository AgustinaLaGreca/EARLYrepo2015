function [getsam, getdsam] = samgetters(D, tau);
% samgetters - functions that get sampled data and their time derivative
%   [getsam, getdsam] = samgetters(D, tau) returns to function handles,
%   getsam and getdsam. When calling these functions with no arguments, as
%   in 
%      sam = getsam()
%      dsam = getdsam(),
%
%   these functions return the concatenated samples in D and their time
%   derivative, respectively. tau is the smooting time used for the
%   derivative.
%
%   This construct can be used to get simple to use handles to the data
%   without causing heavy memory load.
%
%   See also deristats.

getsam = @()local_g(D.ExpID, D.RecID, D.icond);
getdsam = @()local_gd(D.ExpID, D.RecID, D.icond, D.dt_ms, tau);
%getsam = @()smoothen(diff(catsam(readTKABF(D.ExpID, D.RecID, D.icond)))/D.dt_ms, tau, D.dt_ms);

%==================
function [sam, T]=local_g(ExpID, RecID, icond);
[sam, T] = catsam(readTKABF(ExpID, RecID, icond));

function [dsam, T]=local_gd(ExpID, RecID, icond, dt, tau);
[sam, T] = catsam(readTKABF(ExpID, RecID, icond));
dsam = smoothen(diff(sam)/dt, tau, dt);
T = T(1:end-1);






