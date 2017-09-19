function ZColors = getZColors(XYP, plotIndex)
% getZColors       Returns ZColors from XYPlotObject for each plot in
%                plotIndex (defaults to all). 

% Format ZColors 
%	

%% ---------------- CHANGELOG -----------------------
%  Fri Aug 12 2011  Abel   
%   - Initial creation


%% Params 
nrOfPlots = nPlots(XYP);

%% Check Param 
if nargin > 2
    error('getZColors takes two arguments.');
end

if nargin == 1
	plotIndex = 1:nrOfPlots;
end

if plotIndex > nrOfPlots
	error('plot index out of bounds, only %d plots available', nrOfPlots);
end

%% Return ZColors
%ZColors will always be there, in worst case senario its identical to the
%YData
% ZColors = XYP.ZColors(plotIndex);
try 
	ZColors = XYP.ZColors(plotIndex,:);
catch
	ZColors = [];
end
