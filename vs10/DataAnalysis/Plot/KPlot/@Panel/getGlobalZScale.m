function globalZScale = getGlobalZScale(myPanel, plotIndex)
% getZScale      Returns ZScale from Panel for each plot in
%                plotIndex (defaults to all). All ZScales from each XYPlotObject
%                are combined into a single vector. The out value is the
%                sorted, unique values from this vector. 
%
% !ATTENTION, each XYP Object can only have 1 ZScale. If you need more,
% you'll have to add more XYP Objects to the panel.
% FORMAT of the value: Double -> 1 row 


%% ---------------- CHANGELOG -----------------------
%  Fri Aug 12 2011  Abel
%   - Initial creation
%  Tue Nov 29 2011  Abel   
%   - Bugfix + return sorted and unique values 

%% Params
nrOfPlots = nPlots(myPanel);

%% Check Param
if nargin > 2
	error('getZScale takes max two arguments.');
end

if nargin == 1 || strcmpi('all', plotIndex)
	plotIndex = 1:nrOfPlots;
end

if plotIndex > nrOfPlots
	error('plot index out of bounds, only %d plots available', nrOfPlots);
end

%% Return Combined ZScale
globalZScale = {};
zScales = getZScale(myPanel, plotIndex);

% remove empty zScales
zScales = zScales( ~cellfun(@isempty, zScales) );

%return [] if all empty
if isempty(zScales)
	globalZScale = zScales;
	return;

%return single scale if only one scale was found 
elseif size(zScales, 2) == 1 
	globalZScale = zScales{1};
	return;
end

%combine scales 
myTmp = [];
for n=1:length(zScales)
	%Don't include duplicates. The latter are a result of the repmat'ing
	%the scales for each plot (row) within the XYPlotObject
	if n>1 && length(zScales{n-1}) == length(zScales{n}) && all(zScales{n-1} == zScales{n})
		continue;
	end
	myTmp = [ myTmp, zScales{n}];
end
globalZScale = sort(unique(myTmp));

