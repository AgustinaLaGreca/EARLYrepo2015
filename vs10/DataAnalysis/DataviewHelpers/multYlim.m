function MultYlim(axh, YL, X_or_Y);
% MultYlim - impose Y-limits on multiple axes
%    MultYlim(h, [Ymin Ymax]) imposes the same Ylimits on multiple axes with
%    handles h(1), h(2)... .
%
%    MultYlim(h, 'auto') sets the Ylim mode of all the axes to 'auto'.
%
%    MultYlim(h, 'samerange') imposes an identical Y-range to all axes in
%    handle array h. The Y-range is the smallest range that includes all
%    individual Y-ranges when the Ylim mode of h(k) is set to 'auto'.
%
%    MultYlim(h, 'samescale') tolerates differences in offset between the
%    individual Y-ranges, but ensures that the Y-scales are identical. That
%    is, diff(ylim(h(k))) is identical for all k. For log scales, it is the
%    ratio between the limits (rather than the difference) that is made 
%    constant between the individual axes.
%    
%    See also PlotPanes, Ylabels, MultXlim, MultZlim.

if nargin<3, X_or_Y = 'Y'; end

switch upper(X_or_Y),
    case 'X',
        ScaleProp = 'Xscale';
        limfun = @xlim;
    case 'Y',
        ScaleProp = 'Yscale';
        limfun = @ylim;
    case 'X',
        ScaleProp = 'Zscale';
        limfun = @zlim;
end

N = numel(axh);
if isnumeric(YL) || isequal('auto', YL),
    for ii=1:N,
        limfun(axh(ii),YL);
    end
elseif isequal('samerange', YL),
    YL = zeros(N,2);
    % set all to 'auto' and get their limts
    for ii=1:N,
        limfun(axh(ii),'auto');
        YL(ii,1:2) = limfun(axh(ii));
    end
    % now take grand min & max, and impose it to all
    MultYlim(axh,[min(YL(:,1)), max(YL(:,2))],X_or_Y);
elseif isequal('samescale', YL) 
    isLog = isequal({'log'}, unique(get(axh, ScaleProp)));
    if isLog, Yscale='log'; else, Yscale='linear'; end
    YL = zeros(N,2);
    % set all to 'auto' and get their limts
    for ii=1:N,
        set(axh(ii),ScaleProp, Yscale);
        limfun(axh(ii),'auto');
        YL(ii,1:2) = limfun(axh(ii));
    end
    if isLog, YL=log(YL); end
    DY = max(diff(YL,1,2));  % max Y diff
    MY = mean(YL,2); % mean limfun per axes
    YL = [MY-DY/2 MY+DY/2]; % new limits having same DY, but orgina; diff(limfun)
    if isLog, YL=log(YL); end
    for ii=1:N,
        limfun(axh(ii),YL(ii,:));
    end
end







