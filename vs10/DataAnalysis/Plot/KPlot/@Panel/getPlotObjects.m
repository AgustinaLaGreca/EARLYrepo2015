function plotObjects = getPlotObjects(panel, plotIdx)

%% Check params
if nargin == 1
        plotIdx = 1:length(panel.plotObjects);
elseif nargin > 2
	error(['Only two argumens expected.']);
end

plotObjects = panel.plotObjects(plotIdx);

end
