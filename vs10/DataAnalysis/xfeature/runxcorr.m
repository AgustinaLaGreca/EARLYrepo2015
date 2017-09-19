function Y = runxcorr(X, A);
% runxcorr - running normalized crosscorrelation
%   Y = runxcorr(X, A) returns crosscorrelation of arrays X and A,
%   normalized by A and the part of X participating in the xcorr.

X = X(:);
A = A(:);
A = A-mean(A);
A = A/std(A);
A = A(end:-1:1);
M = numel(A);
AvBox = ones(size(A))/M; % boxcar for running averages of X
MeanX = conv(X,AvBox);
%dplot(1, X); xdplot([1 -M/2], MeanX, 'r');
NormX = sqrt((conv(X.^2, AvBox) - MeanX.^2));
Y = conv(X,A)/M;
Y = (Y-MeanX)./NormX;
dplot(1, X); xdplot([1 -M/2], NormX, 'r'); xdplot([1 -M/2], -NormX, 'r');



