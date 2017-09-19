function MultZlim(axh, L);
% MultZlim - impose Z-limits on multiple axes
%    MultZlim(h, [Zmin Zmax]) imposes the same Zlimits on multiple axes with
%    handles h(1), h(2)... .
%
%    MultZlim(h, 'auto') sets the Zlim mode of all the axes to 'auto'.
%
%    MultZlim(h, 'samerange') imposes an identical Z-range to all axes in
%    handle array h. The Z-range is the smallest range that includes all
%    individual Z-ranges when the Zlim mode of h(k) is set to 'auto'.
%
%    MultZlim(h, 'samescale') tolerates differences in offset between the
%    individual Z-ranges, but ensures that the Z-scales are identical. That
%    is, diff(zlim(h(k))) is identical for all k. For log scales, it is the
%    ratio between the limits (rather than the difference) that is made 
%    constant between the individual axes.
%    
%    See also PlotPanes, Ylabels, MultYlim, MultZlim.

multYlim(axh, L, 'Z'); % delegate to multYlim



