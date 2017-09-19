function XYP = setZColors(XYP, ZColors, plotIndex)
%setZColors          Sets XYPlotObject.ZColors in RGB format for each of
%                    the X/Y points.
%
%DESCRIPTION         The colors are set based on: (see doc caxis,
%                    "How Color Axis Scaling Works")
% 			 			index = fix((C-cmin)/(cmax-cmin)*m)+1;
% 						%Clamp values outside the range [1 m]
% 			 			index(index<1) = 1;
% 			 			index(index>m) = m;
%                    Where
%                     cmin = the maximum of XYPlotObject.ZScale
%                     cmin = the minimum of XYPlotObject.ZScale
%                     m = number of colors
%                    The coloring uses either XYPlotObject.ZData or
%                    XYPlotObject.YData (if no ZData was found) to set the
%                    colors.
%INPUT
%                    XYP = XYPlotObject()
%                    ZColors = cell of RGB colors
%                    plotIndex (optional) = idx of XYPlotObject() plots to
%                    set. Defaults to all plots.


%% ---------------- CHANGELOG -----------------------
%  Fri Aug 12 2011  Abel
%   - Initial creation
%  Thu Aug 18 2011  Abel   
%   - Added support for ZNaNColor option
%  Tue Nov 22 2011  Abel   
%   - Partial re-write since ZColors are no longer stored in the Panel()
%     instance. 
%  Tue Nov 29 2011  Abel   
%   - Addapted to new Z-structure within XYPlotObjects
%   - bugfix: NaN color only set for lines not individual points


%% Check Param
nrOfPlots = nPlots(XYP);
nrOfColors = size(ZColors, 1);

if nargin < 2
	error('setZColors takes 2 or 3 arguments.');
elseif nargin == 2
	plotIndex = 1:nrOfPlots;
end

if strcmpi('all', plotIndex);
	plotIndex = 1:nrOfPlots;
end

%% Set Colors within ZScale based on ZData

ZScale = getZScale(XYP);
if isempty(ZScale)
	error('No Zscale set on XYP object, run addZColor on Panel() object first');
end

%ZScale in build as minYval:calculated Y increment:maxYval;
%-> so take first and last as min and max 
minY = ZScale(1);
maxY = ZScale(end);

for cPlot = plotIndex
	%Get ZData from XYPlotObject
 	ZData = cell2mat( torow(getZData(XYP, cPlot)) );

	%Map Z-Value to color index (see doc caxis, "How Color Axis Scaling Works")
	%for each plot
	for n = 1:length(ZData)
		%if NaN, set nan color option for point ZData(n)
		if isnan(ZData(n))
			warning('SGSR:Info', '%s\n', 'Found a NaN value in the Z-Data, the default behaviour is to skip this point',...
			'Optionally you can select a "NaN" color by setting the "ZNaNColor" option',...
			'(see "doc KStructplot" or "doc XYPlotObject")');
			XYP.ZColors{cPlot,n} = cell2mat(get(XYP, 'ZNaNColor'));
			continue;
		end
		%Calculate idx
		index = fix((ZData(n)-minY)/(maxY-minY)*nrOfColors)+1;
		index(index<1) = 1; %Also solves minY == maxY -> index = - Inf
		index(index>nrOfColors) = nrOfColors;
		%set color
		XYP.ZColors{cPlot,n} = ZColors(index, :);
		%Attention!, There will be empty entries in some rows of XYP.ZColors
		%since matlab makes the cell square and fills missing values with blanks
	end
	
end
