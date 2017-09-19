function ZWdemoII(expname, iDataset, datatype);
% ZWdemoII(expname, iDataset, datatype);

[expname, iDataset, datatype] = arginDefaults('expname, iDataset, datatype', 'RG11353', 35, 'phase');

iChan = 2;
iCond = 1;

Norder = 3;             %order of DPs to consider
phaseref = 'response';  %phase referencing of DPs
ME = [];                %transfer function to apply, [] = none

%datatype = 'phase';     %look at phase, see ZW_unraffleDPs.m
bandtoplot = 21:29;     %which equal-Fdp bands to plot, see ZW_plotunraffle.m

%--- import the dataset, if not already ---
try,
    exist(DS);
catch
    DS=[];
end
if isempty(DS), DS = getds(expname); end
tmp = disp(DS);
indx = strfind(tmp,'experiment')+11;
if ~isequal(upper(tmp(indx:end)),upper(expname)),
    DS = getds(expname);
end
    
try,
    exist(ds);
catch,
    ds = [];
end
if isempty(ds), ds = DS(iDataset); end
if ds.ID.iDataset~=iDataset, ds = DS(iDataset); end

    
%--- get the necessary information ---
[df, Aspec, Cspec]=ZWspec(ds, iChan, iCond); %spectra, for plotting
Prim = ZW_getPrim(ds, iChan, iCond, ME); %primaries-info, for plotting
DP = ZW_residualphase(ds, iChan, iCond, ME, Norder, phaseref); %DP information
rDP = ZW_unraffleDPs(DP, datatype); %reorder DPs in equal-Fdp bands

%--- helper for making plotting titles ---
Tstr = [upper(name(ds.ID.Experiment)) ', Rec ' num2str(ds.ID.iDataset)];
    
%--- figure 1: plot of amplitude spectrum ---
iNyquist = round(numel(Aspec)/2); %get the noise floor
noisefreq = linspace(0.01, iNyquist*df,1e3);
noiseampl = ZWnoise(ds, iChan, iCond, ME, noisefreq*1e3);

figure;
dplot(df,Aspec,'k'); %spectrum
xplot(noisefreq, noiseampl,'c'); %plot noisefloor
iKp = DP.wght>0 & ~DP.notUnique & ~DP.notallPrim; %mark significant DPs
iKp = full(iKp);
h(3) = xplot(DP.Fdp(iKp)/1e3, DP.mg(iKp),'gp','markerfacecolor','g','markersize',6);
Ndp = int2str(sum(iKp));
iKp = DP.wght>0 & ~DP.notUnique & DP.notallPrim; %mark significant DPs, but not all prims. signif.
qq = xplot(DP.Fdp(iKp)/1e3, DP.mg(iKp),'bp','markerfacecolor','b','markersize',6);
Lstr = {'sign. prim.' 'non-sig. prim.' 'sign. DP; prim. OK'};
if ~isempty(qq),
    h(4)=qq;
    Lstr = [Lstr {'sign. DP; prim. not OK'}];
end
    
iKp = Prim.wght>0; %mark primaries
h(1) = xplot(Prim.Fprim(iKp)/1e3, Prim.mg(iKp),'rs','markerfacecolor','g','markersize',5);
qq = xplot(Prim.Fprim(~iKp)/1e3, Prim.mg(~iKp),'ro','markerfacecolor','r','markersize',4);
if ~isempty(qq), h(2)=qq; end

xlim([0 df*iNyquist]); %prettify
ylim([-90 -15]);
xlabel('frequency [kHz]');
ylabel('amplitude [dB]');
title([Tstr '; ' Ndp ' significant DPs']);
% legend(h,Lstr,'location','northeast');

%--- figure 2: plot of phase for significant DPs ---
iKp = DP.wght>0 & ~DP.notUnique & ~DP.notallPrim;
F = DP.Fdp(iKp)/1e3;
P = DP.ph(iKp);
R = DP.ph_residual(iKp);
M = DP.ph_main(iKp);
[F, P, R, M] = sortAccord(F, P, R, M, F);

figure;
subplot(2,1,1);
plot(F,P,'gp','markerfacecolor','g','markersize',6);
xplot(F, M,'r');
title([Tstr '; phase (re. ' phaseref ') for significant DPs']);
ylabel('phase [cycle]');
xlim([0 df*iNyquist]);

subplot(2,1,2);
plot(F,R,'gp','markerfacecolor','g','markersize',5);
ylabel('resid. phase [cycle]');
xlabel('frequency [kHz]');
xlim([0 df*iNyquist]);

%--- figure 3: number of significant & unique DPs per equal-Fdp band ---
figure;
subplot(2,1,1);
clear h;
h(1) = plot(rDP.DPinband(:,1),'r-p');
h(2) = xplot(rDP.DPinband(:,2),'b->');
h(3) = xplot(rDP.DPinband(:,3),'k-*');
title([Tstr '; possible ' int2str(sum(rDP.DPinband(:,1))) ' DPs']);
ylabel('count [-]');
legend(h,{'potential' 'unique' 'signif.'});

subplot(2,1,2);
plot(rDP.DPinband(:,3),'k-*');
xlabel('index for equal-Fdp band');
ylabel('count [-]');
title(['signif. # of DPs (N=' int2str(sum(rDP.DPinband(:,3))) ')']);

%--- series of figures with contours for equal Fdp-bands ---
for ii = bandtoplot,
    ZW_plotunraffle(rDP, ii);
end


