function panel = addZColors(panel, plotIndex, doRedraw)
% addZColors adds Z-axies colors to the XYPlotObjects of the given Panel instance, 
% and redraws if asked
%
% panel = addZColors(KPlotObject panel, plotIndex, 'noredraw')
% adds the plot P to KPlotObject instance panel.
%
%      panel : A Panel instance we wish to add a plot object to.
%      plotIndex : (optional, default == 'all') Set to which
%                  XYPlotObjects() the colors are added. 
%        'all' Add colors to all XYObjects in Panel. All data is on single scale. 
%        'single' Add colors to all XYObjects in Panel. Each XYObject has its own scale 
%                 based on its own ZData. 
%        [index] Set colors on certain XYObjects with [index]
% ~ if [plotIndex] -> set colors on certain XY objects, according to
% plotIndex. The scales may be different
%             
%      'noredraw': (optional) If this string is added to the argument list, the
%                  plot is not redrawn. This might speed up things when
%                  constructing large plots.
%
% Returns: updated panel.
%
% Example:
%  >> P = XYPlotObject(0:0.01:2*pi, sin(0:0.01:2*pi));
%  >> panel = Panel;
%  >> panel = addZColors(panel, 'all');

% Created by: Abel Jonckheer 

%% TODO
% 1) Should we decrease inc in calc_ZScale
%   -> No, since we can automatically increase inc by using 'single' scale
% 2) Update documentation
% 3) if ZColor == 'all' -> the ZScale should be equal for each XYPlotObject, regardless if 
%   the XYObjects have different scales.

%% ---------------- CHANGELOG -----------------------
%  Tue Nov 22 2011  Abel   
%   - re-written to include Zscale option 
%  Tue Nov 29 2011  Abel   
%   - Included TODO nr 3

%% Default params 
nPlotObjects = size(panel.plotObjects, 2);
doSetZScale = false();
useSingleScale = false();

%% Check params
switch nargin
	%Take plotIndex as 'all' by default
	case 1
		plotIndex = 1:nPlotObjects;
		doRedraw = 1;
	case 2
		doRedraw = 1;
	case 3
		if isequal('noredraw', doRedraw)
			doRedraw = 0;
		else
			error(['Only two argumens expected, unless an extra ''noredraw'''...
				' is given. Type ''help addZColors'' for more information.']);
		end
		
	otherwise
		error(['This function expects two or three arguments.'...
			' Type ''help addZColors'' for more information.']);
end

%% If seperate Z-axes for each plotObject, loop over objects and return
if strcmpi('all', plotIndex)
	plotIndex = 1:nPlotObjects;
	%all XYObjects share a single ZScale
	useSingleScale = true();
	%save ZColor setting in Panel()
	panel = set(panel, 'ZColor', 'all', 'noredraw');
elseif strcmpi('single', plotIndex)
	for i=1:nPlotObjects
		panel = addZColors(panel, i, 'noredraw');
	end
	%save ZColor setting in Panel()
	panel = set(panel, 'ZColor', 'single', 'noredraw');
	if isequal(1,doRedraw)
		panel = redraw(panel);
	end
	return;
end
%ZScales After this step:
% ~ If 'all' -> all data on single scale 
% ~ If 'single' -> plotIndex = 1 number of which all the ZData should be on
% a single scale 
% ~ if [plotIndex] -> set colors on certain XY objects, according to
% plotIndex. The scales may be different

%% Add ZColors to XYPlotObjects
%Aquire all YData or ZData from each plotObject in Panel()
%ZData is always present in XYObject (when none is given, ydata is copied)
allZData = [];
ZData = getZData(panel, plotIndex);

%Create colormap if not given
%If we already have a ZScale, use this one to
%calculate the colors, otherwise calc them ourselves 
% - get zScale of each XYObject
zScale = getZScale(panel, plotIndex);

% - if 'ZColor' == 'all' we must force a single ZScale for all
%   XYPlotObjects, even if each XYPlotObject was given a different scale.
if useSingleScale && ~cellRowsIdentical_(zScale)
	zScale = {};
	warning('SGSR:Critical', 'Replacing ZScale values of all XYPlotObjects since ''ZColor'' is set to ''all'' (use single ZScale for all XYPlot''s. See help Panel()');
end

% - if all z-scales are empty, calculate our own based on all ZData 
if all(cellfun(@isempty, zScale))
	for n=1:length(ZData)
		allZData = [ allZData, cell2vect_(ZData{n}) ];
	end
	zScale = {calcZscale_(allZData)};
	doSetZScale = true;
	
	%Repmat zScale for each XYObject
	zScale = repmat(zScale, 1, length(plotIndex));
end

%Loop over XYObjects
% - determine ZScale
% - create colormap
% - set colors 
nCount = 1;
for i = plotIndex
	%What if one of the ZScales on the XY objects is not set? 
	% - if a single scale is wanted (plotIndex == 'all') 
	%   => copy other scale of (next) plotIndex
	% - if plotIndex was an idx list, calculate scale based on combined ZData
	% - if length(plotIndex) == 1 (happens when 'single' was set) 
	%   -> calculate based on ZData of single plot
	if isempty(zScale{nCount}) && useSingleScale == true()
		% * find next non-empty scale
		%non-empty idx:
		idx = find(~cellfun(@isempty, zScale));
		%first occurence 
		idx = idx(1);
		%zScale 
		zScale{nCount} = zScale{idx};
		%save scale in XY obj
		doSetZScale = true();
	elseif isempty(zScale{nCount}) && useSingleScale == false()
		%Calc scale based on Zdata of current XY Obj
		zScale{nCount} = calcZscale_( cell2vect_(ZData{nCount}) );
		%save scale in XY obj
		doSetZScale = true();
	end 
	%save scale in XY obj
	if doSetZScale
		panel.plotObjects{i} = setZScale(panel.plotObjects{i}, zScale{nCount});
	end
		
	%create colormap and set colors on XYObject
	zColors = [];
	nrOfColors = length(zScale{nCount});
	
	mapString = sprintf('%s', ...
		['zColors = ', get(panel, 'ZColorMap'), '(', ...
		int2str(nrOfColors), ');']);
	eval(mapString);

	panel.plotObjects{i} = setZColors(panel.plotObjects{i}, zColors);
	
	nCount = nCount +1;
end


%% Redraw if asked
if isequal(1,doRedraw)
	panel = redraw(panel);
end

%% Local functions
function zScale = calcZscale_(YData)
if length(YData) == 1
	zScale = 1;
	return
end

minY = min(YData);
maxY = max(YData);
%- Estimate the increment of the data by the average of the
%  difference between the sorted and unique data points
inc = mean(diff(unique(sort(denan(YData)))));
%safety for inc=0
if inc == 0
	inc = (maxY-minY)/length(YData);
end

%- Check if the inc is pointing in the correct direction
if minY < maxY && inc < 0
	inc = -1*inc;
elseif minY > maxY && inc > 0
	inc = -1*inc;
end
%- Number of colors needed (lin. scale based in average inc)
zScale = minY:inc:maxY;

function vector = cell2vect_(myCell)
myCell = torow(myCell);
vector = [];
for n = 1:length(myCell)
	vector = [ vector, cell2mat( torow(myCell(n))) ];
end

function isIdent = cellRowsIdentical_(inCell)
isIdent = true();
for n=2:length(inCell)
	
	if length(inCell{n-1}) == length(inCell{n}) && all(inCell{n-1} == inCell{n})
		continue;
	else
		isIdent = false();
		return;
	end
	
end
