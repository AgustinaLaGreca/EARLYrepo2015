function YData = getYData(XYP, plotIndex)
% getYData       Returns YData from XYPlotObject for each plot in
%                plotIndex (defaults to all). 


%% ---------------- CHANGELOG -----------------------
%  Fri Aug 12 2011  Abel   
%   - Initial creation


%% Params 
nrOfPlots = nPlots(XYP);

%% Check Param 
if nargin > 2
    error('getYData takes max two arguments.');
end

if nargin == 1
	plotIndex = 1:nrOfPlots;
end

if plotIndex > nrOfPlots
	error('plot index out of bounds, only %d plots available', nrOfPlots);
end

%% Return YData
YData = XYP.YData(plotIndex);
