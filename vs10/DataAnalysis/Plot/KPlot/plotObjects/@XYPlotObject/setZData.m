function XYP = setZData(XYP, ZData, plotIndex)
% setZData       Sets ZData for each plot in plotIndex (defaults to all). 


%% ---------------- CHANGELOG -----------------------
%  Fri Aug 12 2011  Abel   
%   - Initial creation


%% Check Param 
nrOfPlots = nPlots(XYP);

if nargin < 2
    error('setZData takes 2 arguments.');
elseif nargin == 2
 	plotIndex = 1:nrOfPlots;
end

if strcmpi('all', plotIndex);
	plotIndex = 1:nrOfPlots;
end

if plotIndex > nrOfPlots
	error('plot index out of bounds, only %d plots available', nrOfPlots);
end

if ~iscell(ZData)
	ZData = {ZData};
end

%% Set ZData
%Repmat the ZData to match each YData point
nRows = size(ZData, 1);
if nRows < nrOfPlots
	ZData = repmat(ZData, nrOfPlots, 1);
end
%check if nr of Z equals nr of Y for each row (plot)
for cPlot = plotIndex
	XYP.ZData{cPlot, 1} = ZData{cPlot, 1};
	if length(XYP.ZData{cPlot, 1}) ~=  length(XYP.YData{cPlot, 1})
		error('ZData should have the same number of points as the YData');
	end
end
