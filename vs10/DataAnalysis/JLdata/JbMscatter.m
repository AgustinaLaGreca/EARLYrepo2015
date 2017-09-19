function S = JbMscatter(JbM, Fld1, Fld2, varargin);
% JbMscatter - scatter plot of binaural metrics returned by JLbinintMetrix
%   S = JbMscatter(JbM, 'field1', 'field2', ...)
%     ... denote plot args
%    return arg S contains some basic stats.

X = [JbM.(Fld1)];
Y = [JbM.(Fld2)];
S.corr = corr(X(:), Y(:));
S.P_ord1_fit = polyfit(X,Y,1);
S.coeff_linfit = X/Y;

h = xplot(X,Y,varargin{:});
FNS = unique({'iexp' 'icell' 'SPL' 'freq' 'iseqi' 'iseqc' 'iseqb' 'Ui' 'Uc' 'Ub' Fld1 Fld2});
Jr = structpart(JbM, FNS);
clear JbM % save memory; otherwise anonymous fcns below cause memory problems
%dsize(Jr, X, Y)
IDpoints(h, Jr, 1:numel(Jr), @local_label, ...
    'JLbinint', @(Jr,k)JLbinint(Jr(k).Ui, Jr(k).Uc, Jr(k).Ub, 1), ...
    '2D ipsi', @(Jr,k)JLanova(Jr(k).Ui), ...
    '2D contra', @(Jr,k)JLanova(Jr(k).Uc), ...
    '2D bin', @(Jr,k)JLanova(Jr(k).Ub), ...
    '-spont history', @(Jr,k)JLspontHist(Jr(k).Ub), ...
    'raw data', @(Jr,k)local_clampfit(Jr(k).Ub), ...
    '-hist X', @(Jr,k)local_hist(Jr,k, Fld1), ...
    'hist Y', @(Jr,k)local_hist(Jr,k, Fld2) ...
    );

xlabel(Fld1, 'Interpreter', 'none');
ylabel(Fld2, 'Interpreter', 'none');

S.xlim = xlim;
S.ylim = ylim;
S.minLim = min(min(xlim), min(ylim));
S.maxLim = max(max(xlim), max(ylim));

%==========================
function Str = local_label(J, k);
%J, p
p = J(k);
Str = [num2str(p.iexp) '/' num2str(p.icell) '; {' trimspace(num2str([p.iseqi p.iseqc p.iseqb])) '} ' ...
    num2str(p.freq) ' Hz; ' num2str(p.SPL) ' dB' ];

function local_hist(J, k, Fld);
figure;
hist([J.(Fld)], round(numel(J)/25));
xlabel(Fld, 'interpreter', 'none');

function local_clampfit(Uidx);
clampfit(Uidx);
subplot(3,1,1);
R = JLgetRec(Uidx);






