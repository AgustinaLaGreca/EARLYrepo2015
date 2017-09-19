function S = RueCorrPop(L, excludeZerodB);
% RueCorrPop - population analysis of crosscorrelation values
%    S = RueCorrPop(L, excludeZerodB) evaluates the crosscorrelation values computed by
%    RuePopCorr across the cells listed in cell string L. Typically, L is
%    the return vaue of RueList, or a selection thereof. excludeZerodB is a
%    logical deciding whether 0-dB conditions should be excluded. Default
%    is Flase, i.e., include 0-dB conditions.
%
%    The output arg S is a struct array with fields
%          FN1: name of first cell in comparisons
%          FN2: name of second cell in comparisons
%       isPair: logical array indicating whether correlated cells are a true pair.
%       isSame: logical array indicating whether correlated cells are identical.
%
%        Nrep1: # reps for cell 1 
%        Nrep2: # reps for cell 2
%
%       rho_nd: normalized corr averaged over the nondiagonal rep pairs.
%       std_nd: standard deviation of rho_nd
%     ttest_nd: ttest of nondiagonal rho values, testing rho~=0 at 5% conf
%    ttestR_nd: one-tailed ttest of nondiagonal rhos, rho>0 at 5% conf level
%
%     rho_diag: normalized corr averaged over the diagonal rep pairs.
%     std_diag: standard deviation of rho_diag
%   ttest_diag: ttest of diagonal rho values, testing rho~=0 at 5% conf
%  ttestR_diag: one-tailed ttest of diagonal rhos, rho>0 at 5% conf level
%
%      rho_all: normalized corr averaged over the all rep pairs (whether diagonal or not).
%      std_all: standard deviation of rho_all
%    ttest_all: ttest of all rho values, testing rho~=0 at 5% conf
%   ttestR_all: one-tailed ttest of all rhos, rho>0 at 5% conf level
%
%   rho_disatt: disattenuated rho values
%         rho1: mean of nondiagonal rho values of cell 1
%         rho2: mean of nondiagonal rho values of cell 2
%    tt_both_1: conbined ttest of rho1~=0 and rho2~=0, 1% conf level
%    tt_both_5: conbined ttest of rho1~=0 and rho2~=0, 5% conf level
%
%    See also RuePopCorr.

if nargin<2, excludeZerodB=0; end

NN = numel(L);
[rho_nd rho_diag rho_all isPair isSame] = deal([]);
S = [];
for i1=1:NN,
    p1 = RuePopcorr(L{i1}, L{i1} ,excludeZerodB);
    for i2=i1:NN,
        p = RuePopcorr(L{i1}, L{i2}, excludeZerodB);
        p2 = RuePopcorr(L{i2}, L{i2}, excludeZerodB);
        [rho_nd, rho_diag, rho_all] = local_dndfun(@mean, p.rho);
        [std_nd, std_diag, std_all] = local_dndfun(@std, p.rho);
        [ttest_nd, ttest_diag, ttest_all] = local_dndfun(@ttest, p.rho);
        [ttestR_nd, ttestR_diag, ttestR_all] = local_dndfun(@ttest, p.rho,[],[],'right');
        rho1 = local_dndfun(@mean, p1.rho);
        rho2 = local_dndfun(@mean, p2.rho);
        tt1_1 = local_dndfun(@ttest, p1.rho, 0, 0.01,'right');
        tt2_1 = local_dndfun(@ttest, p2.rho, 0, 0.01,'right');
        tt_both_1 = tt1_1&&tt2_1;
        tt1_5 = local_dndfun(@ttest, p1.rho, 0, 0.01,'right');
        tt2_5 = local_dndfun(@ttest, p2.rho, 0, 0.01,'right');
        tt_both_5 = tt1_5&&tt2_5;
        rho_disatt = rho_nd/sqrt(rho1*rho2);

        s1 = structpart(p, {'FN1' 'FN2' 'excludeZerodB' 'isPair' 'isSame' '-' 'Nrep1' 'Nrep2' '-'});
        s2 = CollectInStruct(...
            rho_nd, std_nd, ttest_nd, ttestR_nd, '-', ...
            rho_diag, std_diag, ttest_diag, ttestR_diag, '-', ...
            rho_all, std_all, ttest_all, ttestR_all, '-', ...
            rho_disatt, rho1, rho2, tt_both_1, tt_both_5);
        S = [S structJoin(s1,s2)];
    end
end


% ============
function [Ynondiag, Ydiag, Yall] = local_dndfun(fun, M, varargin);
% apply function fun to diagonal, nondiagonal, and all matrix elements
Npair = numel(M);
Nrep1 = size(M,1);
idiag = 1:(Nrep1+1):Npair; % indices in p.rho and peers selecting diagonal reps
inondiag = setdiff(1:Npair, idiag);
Ydiag = feval(fun, Columnize(M(idiag)), varargin{:});
Ynondiag = feval(fun, Columnize(M(inondiag)), varargin{:});
Yall = feval(fun, Columnize(M), varargin{:});









