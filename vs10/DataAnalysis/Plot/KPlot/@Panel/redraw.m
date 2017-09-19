function panel = redraw(panel)
% REDRAW Redraws the Panel.
%
% redraw(Panel panel)
% This function is mostly used internally, and users only need it when they
% closed the figure and want to see it again.

% Created by: Kevin Spiritus
% Last edited: December 4th, 2006

%% ---------------- CHANGELOG -----------------------
%  Fri Mar 4 2011  Abel   
%   - Added option to keep standard matlab log-scale
%  Tue Apr 19 2011  Abel   
%   - Added legendObject handling 
%  Mon Nov 21 2011  Abel   
%   - Added draw colorbar support
%  Tue Apr 24 2012  Abel   
%   - Added support for x/yaxes false
%%

%% Make sure figure is there
panel = resetAxes(panel);   %sets panel.hdl (= gca handle )

%% Empty the current plot
delete( get(panel.hdl, 'children') );

%% Set the title
axes(getHdl(panel));

title(panel.params.title);
xlabel(panel.params.xlabel);
ylabel(panel.params.ylabel);
if isequal([0 0], panel.params.xlim)
    panel.params.xlim = 'auto';
end
if isequal([0 0], panel.params.ylim)
    panel.params.ylim = 'auto';
end
xlim(panel.params.xlim);
ylim(panel.params.ylim);

%% ticks
mainAxes = gca;

%by Abel: make code compatible with matlab2015b where handles are replaced
%by objects
mainAxes = double(mainAxes);

set(mainAxes, 'TickDir', panel.params.ticksDir, 'Box', panel.params.box);
if ~isempty(panel.params.xTicks)
    set(mainAxes, 'XTick', panel.params.xTicks);
end
if ~isempty(panel.params.xTickLabels)
    set(mainAxes, 'XTickLabel', panel.params.xTickLabels)
end
if ~isempty(panel.params.yTicks)
    set(mainAxes, 'YTick', panel.params.yTicks);
end
if ~isempty(panel.params.yTickLabels)
    set(mainAxes, 'YTickLabel', panel.params.yTickLabels)
end

%% if requested, draw right Y-axes
if isequal('yes', lower(panel.params.rightYAxes))
    rightAxes = axes('Position', get(mainAxes, 'position'), 'TickDir', ...
        panel.params.ticksDir, 'YAxisLocation', 'right', 'XTick', [], ...
        'XTickLabel', '', 'YLim', panel.params.ylim, ...
        'YTick', panel.params.rightYPositions, 'YTickLabel', panel.params.rightYLabels);
    set(get(rightAxes, 'YLabel'), 'String', panel.params.rightYLabel);
    set(mainAxes, 'Color', 'none'); % make sure right axes is visible
end

%% set axes logarithmic if asked
switch(panel.params.logX)
    case 'yes'
        set(mainAxes, 'XScale', 'log');
		%by Abel: Set labels for tics X-as
		if isempty(panel.params.xTickLabels) && strcmpi(panel.params.logXformat, 'xlog125')
			xlog125;
		end
		
    case 'no'
        set(mainAxes, 'XScale', 'linear');
    otherwise
        error('Wrong format for logX');
end

switch(panel.params.logY)
    case 'yes'
        set(mainAxes, 'YScale', 'log');
		%by Abel: Set labels for tics Y-as
		if isempty(panel.params.yTickLabels) && strcmpi(panel.params.logYformat, 'ylog125')
			ylog125;
		end
    case 'no'
        set(mainAxes, 'YScale', 'linear');
    otherwise
        error('Wrong format for logY');
end

%% Disable axes if requested
if ~panel.params.axes
    set(mainAxes, 'Visible', 'off');
end
if ~panel.params.xaxes
	set(mainAxes, 'XTick', []);
end
if ~panel.params.yaxes
	set(mainAxes, 'YTick', []);
end


%% Reverse axes if requested
if isequal(panel.params.reverseY, 'yes')
    set(mainAxes, 'YDir', 'reverse');
end
if isequal(panel.params.reverseX, 'yes')
    set(mainAxes, 'XDir', 'reverse');
end

%% Set Tag if requested
if ~isempty(panel.params.Tag)
	set(mainAxes, 'Tag', panel.params.Tag);
end

%% Then, plot all plot objects
nPlotObjects = size(panel.plotObjects, 2);
% remember the colors, so all plots can get their own colors
startColor = 1;
for i = 1:nPlotObjects
    panel.plotObjects{i} = redraw(mainAxes, panel.plotObjects{i}, startColor);
    startColor = startColor + nPlots(panel.plotObjects{i});
end

nTextBoxes = size(panel.textObjects, 2);
for i = 1:nTextBoxes
    panel.textObjects{i} = redraw(mainAxes, panel.textObjects{i});
end

%by Abel: draw colorbox if requested
if get(panel, 'addZColorbar') && strcmpi(get(panel, 'ZColor'), 'all'); 
 	%get ZScale & Colors
	zScale = getGlobalZScale(panel);
	minZ = min(zScale);
	maxZ = max(zScale);
	set(mainAxes, 'CLim', [minZ, maxZ]);
	
	%calc colors for global Z scale
	zColors = [];
	nrOfColors = length(zScale);
	
	mapString = sprintf('%s', ...
		['zColors = ', get(panel, 'ZColorMap'), '(', ...
		int2str(nrOfColors), ');']);
	eval(mapString);

	%draw bar
	%Strange matlab behaviour if nr of colors > 200 
	%-> in this case we need to create the colorbar ourselves
	if length(zScale) > 200
		%let matlab build a colorbar axes and remove its colors
		cbh = colorbar('peer', mainAxes);
		delete(get(cbh, 'Children'));
		
		%start drawing patches 
		%-> Calculate color increments
		%- get nr of yticks in colorbar axes
		yticks = get(cbh, 'YTick');
		%- we want 10 sub ticks between the axes increments
		totNrOfSubs = 10 * length(yticks);
		%-> Create coordinates for patch
		%- Coordinates Axes: x = XLim, y = [min(ZScale) max(ZScale)]
		yCoords = linspace(minZ, maxZ, totNrOfSubs);
		xCoords = get(cbh, 'XLim');
		%Calculate color idx
		colorIdx =[];
		for n = 1:length(yCoords)
			%Calculate idx
			colorIdx(n) = fix((yCoords(n)-minZ)/(maxZ-minZ)*nrOfColors)+1;
		end
		
		%set colobar axes active
		axes(cbh);
		%loop over color patches 
		for nRow = 1:length(yCoords)-1
			patch([xCoords(1) xCoords(1) xCoords(2) xCoords(2)],...
				[yCoords(nRow), yCoords(nRow+1), yCoords(nRow+1), yCoords(nRow)],...
				'w', 'FaceColor', zColors(colorIdx(nRow),:), 'EdgeColor', 'none');
		end
		%set active axes back to mainAxes 
		axes(mainAxes);
	else
		colormap(zColors);
		colorbar('peer',mainAxes);
	end
	
elseif get(panel, 'addZColorbar');
	warning('SGSR:Critical', 'Can''t plot a colorbar since the ZScale is not shared over all elements of the Panel. Set the ZColor option of Panel() to ''all''');
end

%by Abel: add legend objects (must be last??->NOT)
nLegends = size(panel.legendObjects, 2);
for i = 1:nLegends
    panel.legendObjects{i} = redraw(panel.legendObjects{i}, panel);
end

%% Add plots to the UserData field of the page
if ~panel.isRegistered
    pages = findobj('Tag','defaultpage');
    if ~isempty(pages)
        page = pages(1);
        currentChilds = length(get(page, 'UserData'));
        if length(page) > 1
            for n = 2:length(pages)
                numChilds = length(get(pages(n), 'UserData'));
                if numChilds < currentChilds
                    page = pages(n);
                end
            end
        end
        userdata = get(page,'UserData');
        panel.isRegistered = true;
        if iscell(userdata)
            userdata{length(userdata)+1} = panel;
        else
            userdata{1} = panel;
        end
        set(page,'UserData',userdata);
    end
end
