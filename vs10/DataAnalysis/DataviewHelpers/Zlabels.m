function Zlabels(axh, Txt, varargin);
% Zlabels - put identical Z-labels to multiple axes.
%    Zlabels([h1 h2 ..], Str) puts label Str besides the Z-axis of axes h1,
%    h2, .. . 
%
%    Zlabels([h1 h2 ..], Str, Prop1, Val1, ...) also passes property/value
%    pairs to ZLABEL. 
%
%    
%    See also ZLABEL, PlotPanes, Xlabels.

for h=axh,
    zlabel(h,Txt,varargin{:});
end








