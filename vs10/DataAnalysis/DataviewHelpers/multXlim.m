function MultXlim(axh, L);
% MultXlim - impose X-limits on multiple axes
%    MultXlim(h, [Xmin Xmax]) imposes the same Xlimits on multiple axes with
%    handles h(1), h(2)... .
%
%    MultXlim(h, 'auto') sets the Xlim mode of all the axes to 'auto'.
%
%    MultXlim(h, 'samerange') imposes an identical X-range to all axes in
%    handle array h. The X-range is the smallest range that includes all
%    individual X-ranges when the Xlim mode of h(k) is set to 'auto'.
%
%    MultXlim(h, 'samescale') tolerates differences in offset between the
%    individual X-ranges, but ensures that the X-scales are identical. That
%    is, diff(xlim(h(k))) is identical for all k. For log scales, it is the
%    ratio between the limits (rather than the difference) that is made 
%    constant between the individual axes.
%    
%    See also PlotPanes, Xlabels, MultYlim, MultZlim.

multYlim(axh, L, 'X'); % delegate to multYlim


