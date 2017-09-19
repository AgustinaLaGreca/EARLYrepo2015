function Xlabels(axh, Txt, varargin);
% Xlabels - put identical X-labels to multiple axes.
%    Xlabels([h1 h2 ..], Str) puts label Str besides the X-axis of axes h1,
%    h2, .. . The axes handle array is typically the Bh output arg of
%    plotPlanes.
%
%    Xlabels([h1 h2 ..], Str, Prop1, Val1, ...) also passes property/value
%    pairs to XLABEL.
%
%    
%    See also XLABEL, PlotPanes, Ylabels.

for h=axh,
    xlabel(h,Txt,varargin{:});
end








