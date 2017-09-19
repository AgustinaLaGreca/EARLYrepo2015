function YData = getYData(myPanel, plotIndex)
% getYData       Returns YData from Panel for each plot in
%                plotIndex (defaults to all). 


%% ---------------- CHANGELOG -----------------------
%  Fri Aug 12 2011  Abel   
%   - Initial creation


%% Params 
nrOfPlots = nPlots(myPanel);

%% Check Param 
if nargin > 2
    error('getYData takes two arguments.');
end

if nargin == 1
	plotIndex = 1:nrOfPlots;
end

if plotIndex > nrOfPlots
	error('plot index out of bounds, only %d plots available', nrOfPlots);
end

%% Return YData
YData = {};
for n = plotIndex
	YData{end+1} = getYData(myPanel.plotObjects{n});
end
