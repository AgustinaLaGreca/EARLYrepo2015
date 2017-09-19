function DX = derivative(X, M);
% Derivative - estimated time-derivative of sampled signal
%   D=Derivative(X, M) estimates the time derivative D of array X, whereby
%   X is considered as a sampled time signal. M is the length of each the
%   tails used for convolution. Default M = 10. D has one sample less than
%   X.
%
%   See also Smoothen, Diff.

if nargin<2, M=10; end

ii=-M:M; ii(M+1) = [];
Win = ((-1).^(ii+M+1)./ii).';

if isvector(X) && (size(X,1)==1), % row vector -> temp turn into col vector
    X = X(:); 
    isRow=1;
else,
    isRow=0;
end

DX = convmult(X, Win);
DX = DX(1+M:end-M);

if isRow, DX = DX.'; end




