function ZData = getZData(XYP, plotIndex)
% getZData       Returns ZData from XYPlotObject for each plot in
%                plotIndex (defaults to all). 


%% ---------------- CHANGELOG -----------------------
%  Fri Aug 12 2011  Abel   
%   - Initial creation


%% Params 
nrOfPlots = nPlots(XYP);

%% Check Param 
if nargin > 2
    error('getZData takes two arguments.');
end

if nargin == 1
	plotIndex = 1:nrOfPlots;
end

if plotIndex > nrOfPlots
	error('plot index out of bounds, only %d plots available', nrOfPlots);
end

%% Return ZData
%ZData will always be there, in worst case senario it's identical to the YData
ZData = XYP.ZData(plotIndex);
