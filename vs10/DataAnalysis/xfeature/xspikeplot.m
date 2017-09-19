function xspikeplot(X);
% xspikeplot - plot spikes extracted by xspikeX.
%   xspikeplot(X);

[Ev, Tev] = getevent(X);

figure; set(gcf,'units', 'normalized', 'position', [0.612 0.129 0.305 0.177])
hist(X.ymax(X.catIdx>0),100);

Template = [];
figure; clf; set(gcf,'units', 'normalized', 'position', [0.46641 0.38574 0.46953 0.54004])
LegStr = {};
for icat=1:X.Ncat,
    iok = X.catIdx==icat;
    if any(iok),
        xplot(Tev, Ev(:,iok), lico(icat,gca));
        %xplot(Tev, Ev(:,iok));
    end
    Template = [Template, median(Ev(:,iok), 2)];
    LegStr = [LegStr, ['N=' num2str(sum(iok))]];
end

figure; clf; set(gcf,'units', 'normalized', 'position', [0.0015625 0.3925 0.4375 0.525]);
%dsize(Tev, X.tshift)
plot(Tev, Template,'linewidth',2);
legend(LegStr);
grid on
xplot(xlim, [1 1]*X.AlignThr, 'k:');











