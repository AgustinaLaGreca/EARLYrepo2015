function myPanel = setZScale(myPanel, zScale, plotIndex)
%setZScale       Sets ZScale for each plot in plotIndex (defaults to all). 
%   setZscale(Panel) is able to set separate ZScales for each XYObject. However, it can not set
%   different ZScales on each plot (row) within an XYObject.


%% ---------------- CHANGELOG -----------------------
%  Fri Aug 12 2011  Abel   
%   - Initial creation 

%% Check Param 
nrOfPlots = nPlots(myPanel);

if nargin < 2
    error('setZColors takes 2 arguments.');
elseif nargin == 2
 	plotIndex = 1:nrOfPlots;
end

if strcmpi('all', plotIndex);
	plotIndex = 1:nrOfPlots;
end

if plotIndex > nrOfPlots
	error('plot index out of bounds, only %d plots available', nrOfPlots);
end

if ~iscell(zScale)
	zScale = {zScale};
end

%% Set ZScale
%Repmat the zScale to match the nr of XYPlotObjects if needed. 
nRows = size(zScale, 1);
if nRows < nrOfPlots
	zScale = repmat(zScale, nrOfPlots, 1);
end

%Save ZScale in the XYPlotObject
for cPlot = plotIndex
	myPanel.plotObjects{cPlot} = setZScale(myPanel.plotObjects{cPlot}, zScale(cPlot))
end
