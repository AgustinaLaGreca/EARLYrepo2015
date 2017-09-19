function [T, R, D] = MultiMedian(M, I, Nsmooth);
% MultiMedian - medians across various columns of a matrix.
%   TM = MultiMedian(M, I), where M is a matrix, returns the templates TM
%   of subsets of columns specified by cell array I:
%
%       TM(:,k) = median(M(:, I{k}),2)
%
%   Nans are returned for empty I{k}.
%
%   [TM, R, D] = MultiMedian(M, I, Nsmooth) also returns the residues, ie, 
%   M with for each column the corresponding smoothed median subtracted.
%   D(k) is the normalized distance (see NormDist) of M(:,k) with the 
%   median of its class. Nsmooth is the order of the hamming window used 
%   for smoothing (see Smoothen). Default Nsmooth = 10;
%
%   Note:
%   I is typically the output of Categorize.
%
%   See also Categorize, CatPlot, MEDIAN.

if nargin<3, Nsmooth = 10; end

R = M;
D = nan(1,size(M,2));
hasCol = any(cellfun(@(x)size(x,1), I)>1); % true is I contains col arrays
if hasCol,
    Tgrand = median(M(:,cat(1,I{:})),2);
else,
    Tgrand = median(M(:,[I{:}]),2);
end
for k=1:numel(I),
    mm = median( M(:,I{k}), 2);
    T(:,k) = mm;
    D(I{k}) = NormDist(M(:,I{k}), mm);
    N = numel(I{k}); 
    Tsc = Tgrand+smoothen(mm-Tgrand,Nsmooth);
    %f8; plot(Tsc); pause;
    R(: , I{k}) = R(: , I{k})- repmat(Tsc, 1, N); 
end
%R = R-samesize(median(R,2), R);

% if nargout<3, return; end
% 
% % ===advanced residues, obtained by a linear combination of templates
% Nt = size(T,1); % # samples in each wave = column
% TT = smoothen(T,Nsmooth);
% [CC, SS] = princomp(TT); % coefs, scores
% SS = SS./samesize(sqrt(sum(SS.^2)), SS); % normalize all columns
% SS = [SS, ones(Nt,1)/Nt]; % allow a DC comp
% %dsize(SS, M);
% %SS = SS(:,1:3);
% Mapprx = SS*SS.'*M; % approximation of waves in M using lin comb of columns in SS
% avR = M-Mapprx;
% avR = avR-samesize(median(avR,2), avR);







