function [I, Dist] = locateSegment(Seg, X);
% locateSegment - locate segment within array
%    I = locateSegment(Seg, X) finds the position I in X for which array
%    Seg and X(I+(1:N)) are optimally close in a Euclidian sense, where
%    N is length(Seg). If X is a matrix, I is a row vector containing the
%    locations if Seg with each of the columns of X.
%
%    [I, D] = locateSegment(Seg, X) also returns the euclidean distance D
%    between the Seg and the nearest segment within X.
%
%    See also maxcorr, FIND.

if isvector(X), X = X(:); end % force into col vector
Lseg = numel(Seg);
NX = size(X,1);

% compute <Seg*Xseg> term, with Xseg a shifting segment of X
CM = convmtx(Seg(:)',NX-Lseg+1);
SX = CM*X;
% compute <Xseg*Sseg> term
CM = convmtx(ones(1,Lseg),NX-Lseg+1);
XX = CM*X.^2;
% segment^2
SS = sum(Seg.^2);
% distance^2 (Seg-Xseg)^2 = SS + XX - 2SX
D2 = XX + SS - 2*SX;
[D2, I] = min(D2);
Dist = sqrt(D2);

