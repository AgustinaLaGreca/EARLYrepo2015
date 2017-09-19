%===========JS grapes===========
% preliminary analyses of mask data


% SACs
DS = read(dataset, 'RG11302',0);
ds = DS(25); 
SPT = spiketimes(ds,1,'no-unwarp'); 
clear XC;
for ii=1:ds.Stim.Presentation.Ncond, 
    [XC(ii,:) tau] = sptcorr(SPT(ii,:), 'nodiag', 7, 0.1);  
end
plot(tau, XC);
xplot(tau, mean(XC,1), 'k', 'linewidth', 3);

% correlograms within spike trains
ds = DS(25); 
XC = 0 ;
SPT = spiketimes(ds,1,'no-unwarp'); 
for ii=1:ds.Stim.Presentation.Ncond, 
    for jj=1:ds.Stim.Presentation.Nrep, 
        [xc tau] = sptcorr(SPT(ii,jj), SPT(ii,jj), 7, 0.1); 
        XC = XC + xc; 
    end; 
end
plot(tau, XC ,'-*')

% inter-spike intervals
ds = DS(25); 
NN = 0 ;
SPT = spiketimes(ds,1,'no-unwarp', 0.25); 
for ii=1:ds.Stim.Presentation.Ncond, 
    for jj=1:ds.Stim.Presentation.Nrep, 
        isi = diff(SPT{ii,jj});
        [N,DT] = hist(isi, 0:0.2:50);
        NN = NN + N;
    end; 
end
bar(DT, NN); xlim([0 49])

% raw data
SPT = spiketimes(ds,1,'no-unwarp', 0); 
[NW, dt t0] = anadata(ds,1,1,1);
subplot(2,1,1); dplot([dt t0], NW); 
fenceplot(SPT{1,1},ylim, 'r')
subplot(2,1,2); dplot(dt, diff(smoothen(NW,0.1,dt))/dt)
set(gcf,'units', 'normalized', 'position', [0.0227 0.489 0.974 0.41])
TracePlotInterface(gcf)




%==================================