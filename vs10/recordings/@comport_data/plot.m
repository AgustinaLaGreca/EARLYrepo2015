function plot(CP, varargin)
% comport_data/plot - PLOT for comport_data objects.
%
%   See comport_data.

set(gcf,'units', 'normalized', 'position', [0.209 0.549 0.753 0.353])
plot(CP.Data.Times/1e3, 100*CP.Data.Samples/512, 'k-');
xplot(CP.Data.Times/1e3, 100*CP.Data.Samples/512, 'dk', 'markerfacecolor', 'k');
xplot(CP.Data.Times/1e3, 100*CP.Data.Samples/512, 'ow', 'markersize', 2, 'markerfacecolor', 'y');
xlabel('Time since start of recording (s)', 'fontsize', 12);
ylabel('Reflectance (%)', 'fontsize', 12);
ylim([0 100]);


