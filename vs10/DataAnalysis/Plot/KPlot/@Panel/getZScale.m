function ZScale = getZScale(myPanel, plotIndex)
% getZScale      Returns ZScale from Panel for each plot in
%                plotIndex (defaults to all).
% !ATTENTION, each XYP Object can only have 1 ZScale. If you need more,
% you'll have to add more XYP Objects to the panel.
% FORMAT of the value: Double -> 1 row 


%% ---------------- CHANGELOG -----------------------
%  Fri Aug 12 2011  Abel
%   - Initial creation


%% Params
nrOfPlots = nPlots(myPanel);

%% Check Param
if nargin > 2
	error('getZScale takes max two arguments.');
end

if nargin == 1
	plotIndex = 1:nrOfPlots;
end

if plotIndex > nrOfPlots
	error('plot index out of bounds, only %d plots available', nrOfPlots);
end

%% Return ZScale
ZScale = {};
for n = plotIndex
	ZScale{end+1} = getZScale(myPanel.plotObjects{n});
end
