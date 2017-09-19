function JLinterpeak(Uidx,ipeak);
% JLinterpeak(Uidx, ipeak);

ipeak = arginDefaults('ipeak',1:2);

[T, W, R] = JLviewrec(Uidx);


RI = R.RecMatIpsi;
RC = R.RecMatContra;
dti = R.dt_Ipsi;
dtc = R.dt_Contra;
while size(RI,1)>200, 
    RI = resample(RI,1,2);
    dti = 2*dti;
end
while size(RC,1)>200, 
    RC = resample(RC,1,2);
    dtc = 2*dtc;
end

% mDCi = nanmean(RI,1); %slow mean per cycle
% mDCc = nanmean(RC,1); %slow mean per cycle
% RI = bsxfun(@minus, RI, mDCi);
% RC = bsxfun(@minus, RC, mDCc);

mWi = nanmean(RI,2); % mean across cycles
mWc = nanmean(RC,2); % mean across cycles


set(figure,'units', 'normalized', 'position', [0.0891 0.139 0.704 0.762])% IPSI
subplot(2,2,1);
dplot(dti, RI);
xdplot(dti, mWi, 'w', 'linewidth',3);
xdplot(dti, mWi, 'k', 'linewidth',2);
[Tpeak, Y, isort] = peakfinder(dti, mWi, 0.1);
fenceplot(Tpeak,ylim, 'k');
if numel(isort)==1,
    tp = Tpeak(isort);
else,
    tp = Tpeak(isort(ipeak));
end
xlim([0 dti*numel(mWi)]);
fenceplot(tp,ylim, 'b', 'linewidth',3);
t2isam = @(t)1+round(t/dti);
% CONTRA
subplot(2,2,2);
dplot(dtc, RC);
xdplot(dtc, mWc, 'w', 'linewidth',3);
xdplot(dtc, mWc, 'k', 'linewidth',2);
[Tpeak, Y, isort] = peakfinder(dtc, mWc, 0.1);
fenceplot(Tpeak,ylim, 'k');
if numel(isort)==1,
    tp = Tpeak(isort);
else,
    tp = Tpeak(isort(ipeak));
end
xlim([0 dtc*numel(mWc)]);
fenceplot(tp,ylim, 'r', 'linewidth',3);
t2isam = @(t)1+round(t/dtc);
title(JLtitle(Uidx))

subplot(2,2,3);
RI = [RI(:,1:2:end-1); RI(:,2:2:end) ]; RI(:,any(isnan(RI))) = [];
%ccc = cov(RI.');
ccc = corr(RI.');
tau = Xaxis(RI(:,1), dti);
h = pcolor(tau, tau, ccc); set(h, 'edgecolor','none'); colorbar
caxis([-1 1]);
xlabel('Time in cycle (ms)')
ylabel('Time in cycle (ms)')
title('IPSI')


subplot(2,2,4);
RC = [RC(:,1:2:end-1); RC(:,2:2:end) ]; RC(:,any(isnan(RC))) = [];
% ccc = cov(RC.');
ccc = corr(RC.');
tau = Xaxis(RC(:,1), dtc);
h = pcolor(tau, tau, ccc); set(h, 'edgecolor','none'); colorbar
caxis([-1 1]);
xlabel('Time in cycle (ms)')
ylabel('Time in cycle (ms)')
title('CONTRA')



