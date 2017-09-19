function CM  =RueCmap;
% RueCmap - colormap for Rue's correlation plots
%   RueCmap runs from red via white to black;
%
%   See also Rueceplot.

G = gray;
POS = G(end:-2:2,:);
NEG = G(1:2:end-1,:);
NEG(:,1) = 1;
CM = [NEG; POS];





