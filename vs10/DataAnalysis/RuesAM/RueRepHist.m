function RueRepHist(P);
% RueRepHist - histogram of across-rep correlations
%    RueRepHist(P)

idiag = 1:(1+size(P.rho,1)):numel(P.rho);  
Pnd = P.rho; Pnd(idiag)=[]; 
Pd = P.rho(idiag); 
hist(Pnd); 
hh = findobj(gca,'type', 'patch'); 
hold on; 
hist(Pd); 
set(hh, 'facecol', 'r');
xlim([-0.5 1]);
title([P.FN1 ' / ' P.FN2]);





