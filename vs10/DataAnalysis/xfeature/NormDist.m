function D = NormDist(X,Y);
% NormDist - nonsymmetric normalized distance
%   D=NormDist(X, Y) is the Euclidean distance between arrays X and Y,
%   normalized by the standard deviation of Y, i.e.,
%
%       D = mean((X-Y).^2)/std(Y);
%
%   If X is a matrix, D(k) equals NormDist(X(:,k),y).
%
%   See also MultiMedian.

if isempty(X) || isempty(Y),
    D = nan(1,size(X,2));
else,
    D = mean((X-SameSize(Y,X)).^2)./std(Y);
end


