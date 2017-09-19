function S = JLtonePCA(ExpID, RecID, icond, Nmax, plotArg);
% JLtonePCA - PCA of tone responses (no beats)

if nargin<2, RecID=[]; end
if nargin<3, icond=[]; end
if nargin<4, Nmax=10; end
if nargin<5, plotArg=''; end

Tign = 20; % first Tign ms ignored

if numel(icond)>1, % combine multiple plots in one graph using recursion
    for ii=1:numel(icond),
        S(ii) = JLtone2(ExpID, RecID, icond(ii), Nmax, ploco(ii));
    end
    return;
end

%---------single condition from here--------------

[D,DS]=readTKABF(ExpID, RecID, icond); 
W = cat(2, D.AD(2:end,1).samples); 
rho = shufcorr(W);
SPL = unique(DS.SPL);
Fcar = DS.Fcar;
Freq = round(Fcar(icond,1));
dt = D.dt_ms;
[Co, We, Ap] = mypca(W,Nmax); 
dplot(dt, Co);


set(gcf,'units', 'normalized', 'position', [0.0344 0.516 0.902 0.384]);
TracePlotInterface(gcf); ylim auto; ylim(ylim);
title([DS.title '  '  num2str(Freq) ' Hz   ' trimspace(num2str(SPL))  ' dB  ; rho = ' num2str(rho,3)]);


