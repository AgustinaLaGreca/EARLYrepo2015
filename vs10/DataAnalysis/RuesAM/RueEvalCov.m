function S = RueEvalCov(S);
% RueSplitCov - evaluate covariance of RueXXXcorr output.
%   T = RueSplitCov(S), where S is the output of RueXXXcorr, returns a
%   struct having fields
%   ...
%   See also ReadXXXcorr.

% collect relavant data from S and rearrange it


Scov = cat(3,S.urho); % covariances of single stim in 3-dim array. 
% Scov(k,i1,i2) is covariance across cells S.FN1 and S.FN2 of stimulus # k,
% rep # i1 of cell 1 and rep # i2 of cell 2.
Scov = mean(Scov, 3); % this is the within-stim contribution to the covariance
% Now Scov(i1,i2) is the total within-sti mcovariance between reps i1 and
% i2 of cell 1 & 2, resp.
hist(Scov(:))


