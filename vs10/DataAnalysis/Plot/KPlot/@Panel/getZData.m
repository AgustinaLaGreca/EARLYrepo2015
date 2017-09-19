function ZData = getZData(myPanel, plotIndex)
% getZData       Returns ZData from Panel for each plot in
%                plotIndex (defaults to all). 


%% ---------------- CHANGELOG -----------------------
%  Fri Aug 12 2011  Abel   
%   - Initial creation


%% Params 
nrOfPlots = nPlots(myPanel);

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
ZData = {};
for n = plotIndex
	ZData{end+1} = getZData(myPanel.plotObjects{n});
end
