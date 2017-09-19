function [TMPL, T, Nev, Ie] = Mtemplate(M, FN, Bounds, doPlot, FNleg, FNunit);
% Mtemplate - event templates based on spike metrics
%   Mtemplate(M, 'Foo', B) produces templates of events in struct M, the output
%   of spikeMetric, based on spike metric Foo. The templates are
%   across-event medians of groups of events whose Foo values are between
%   the bounds specified in B.
%      Template # 1:  B(1)< M.Foo <=B(2)
%      Template # 2:  B(2)< M.Foo <=B(3), etc.
%
%   [TMPL, T, N, Ie] = Mtemplate(...) also returns time axis T, the number
%   of events underlying each template in TMPL, and their indices Ie.
%
%   Mtemplate(M, 'Foo', B, 1, 'foo', 'AU') also plots the templates, using the
%   name 'foo' for the quantitity in M.(FN) and units AU.
%
%   Mtemplate(M, 'Foo', B, 2, ...) also plots the individual events
%   contained in each category.
%
%   See also SpikeMetrics, getEvent.
if nargin<3,
    Bounds = [-inf inf];
end
if nargin<4,
    doPlot = 0;
end
if nargin<5,
    FNleg = FN;
end
if nargin<6,
    FNunit = '';
end


MT = M.(FN); % the values on which the classification is based
Ncat = numel(Bounds)-1; % # categories defined by boundaries
Nt = numel(getevent(M.X,1)); % # time points per event
TMPL = nan(Nt, Ncat);
Ie = {};
for ic=Ncat:-1:1,
    b0 = Bounds(ic);
    b1 = Bounds(ic+1);
    ihit = betwixt(MT, b0, b1);
    ievent = M.Ie(ihit);
    [Ev, T] = getevent(M.X, ievent);
    tmpl = median(Ev,2);
    if ~isempty(tmpl), TMPL(:,ic) = tmpl; end
    Nev(1,ic) = numel(ievent);
    LegStr{ic} = [num2str(b0) '<' FNleg '<' num2str(b1) ' ' FNunit ' (N=' num2str(Nev(1,ic)) ')'];
    Ie = [Ie, ievent(:)'];
end
if abs(doPlot)>=1,
    figure;
    plot(T,TMPL,'linewidth',2); 
    h=legend(LegStr, 'Location', 'NorthWest'); set (h, 'fontsize',8)
    title([M.ExpID '/' M.RecID '/' num2str(M.icond) ' ---- ' num2str(M.Xval) ' ' M.Xunit]);
    set(gcf,'units', 'normalized', 'position', [0.28 0.336 0.668 0.535]);
    xlabel('Time (ms)', 'fontsize', 12)
    ylabel('V_{ex} (mV)', 'fontsize', 12)
    ff
end
if abs(doPlot)>=2,
    hfig1 = gcf;
    figure; set(gcf,'units', 'normalized', 'position', [0.383 0.0313 0.668 0.535])
    if doPlot<0, Iplot = Ncat-1:-1:1;
    else, Iplot = 1:Ncat-1;
    end
    for iplot=Iplot,
        [Ev, T] = getevent(M.X, Ie{iplot});
        xplot(T, Ev, lico(Ncat-iplot));
    end 
    h=legend(LegStr, 'Location', 'NorthWest'); set (h, 'fontsize',8)
    title([M.ExpID '/' M.RecID '/' num2str(M.icond) ' ---- ' num2str(M.Xval) ' ' M.Xunit]);
    xlabel('Time (ms)', 'fontsize', 12)
    ylabel('V_{ex} (mV)', 'fontsize', 12)
    figure(hfig1); set(gcf,'units', 'normalized', 'position', [0.0367 0.388 0.668 0.535]);
end






