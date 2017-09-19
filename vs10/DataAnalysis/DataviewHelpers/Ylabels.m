function Ylabels(axh, Txt, varargin);
% Ylabels - put identical Y-labels to multiple axes.
%    Ylabels([h1 h2 ..], Str) puts label Str besides the Y-axis of axes h1,
%    h2, .. . The axes handle array is typically the Lh output arg of
%    plotPlanes.
%
%    Xlabels([h1 h2 ..], Str, Prop1, Val1, ...) also passes property/value
%    pairs to YLABEL. 
%
%    
%    See also YLABEL, PlotPanes, Xlabels.

for h=axh,
    ylabel(h,Txt,varargin{:});
end








