function ZColor = getZColors(myPanel, plotIndex)
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
ZColor = {};
for n = plotIndex
	ZColor{end+1} = getZColors(myPanel.plotObjects{n});
end
