function S = anzwuis(Fn, doPlot);
%  anzwuis - zwuis analyse of pilot AN data of Philip
%   S = anzwuis('myfile', doPlot);

load([Fn '_rdam']); % converted data in dataset d
Z = zwuisana(d.Stim, digichan(d,1), 0,0); % 0,0,prep = no baseline, all zwuis cycles, only prepare

rc = rayleigh_rmin(sum(Z.R),0.001);
rc2 = rayleigh_rmin(sum(Z.R),0.0001);
Alpha = rayleigh_pval(abs(Z.Zsp), sum(Z.R));
VS = abs(Z.Zsp);
f2i = @(f)1+round(f/Z.df); % Hz => idx
i2f = @(i)Z.df*(i-1);






