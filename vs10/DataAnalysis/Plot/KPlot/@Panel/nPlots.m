function nPlots = nPlots(panel)
% NPLOTS returns the amount XYPlotObjects in the given panel object
%
% nPlots = nPlots(panel)

% Created by: Kevin Spiritus
% Last edited: December 4th, 2006

nPlots = size(panel.plotObjects, 2);