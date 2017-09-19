function XYPlotObject = redraw(hdl, XYPlotObject, startColor)
% REDRAW Redraws the given plot on the given handle.
%
% XYPlotObject P = Redraw(XYPlotObject P, dbl hdl[, int startColor])
% The plot P is being plotted on the given figure or axes.
%
%           P: The XYPlotObject which is being plotted.
%         hdl: The handle where XYPlotObject is being plotted (typically a figure
%              or axes).
%  startColor: Determines at which position in the ColorOrder array (see
%              axes properties in Matlab Help) the colors of the plot
%              start. This helps making sure different plots have different
%              colors.

% Created by: Kevin Spiritus


%% ---------------- CHANGELOG -----------------------
%  Fri Aug 12 2011  Abel   
%  - Added Z-Axes coloring of data
%  - General code cleanup and documentation
%  Tue Nov 29 2011  Abel   
%  - Bugfix: ZColor not drawn when obj was first drawn without zcolors and
%  then redrawn with zcolors
%  Thu Dec 1 2011  Abel   
%  - Bugfix: Crashes when XYP has no data


%% parameters
if isequal(2, nargin)
	startColor = 1;
end

if nargin > 3
	error('Too many arguments.');
end

%by Abel: check for empty X/Y-values and return if empty 
% - if empty, ML parameters will be empty cell due to repmat while creating
%   the XYP object. This will result in index error since XYPlotObject.params.ML.Color{1}
%   is being checked.
if isempty(XYPlotObject.XData)
	warning('SGSR:Info', 'Nothing to plot in this XYPlotObject: no X-values. =>>> Return to base');
	return;
end
	


%by Abel: Extra params for ZColor
doZColorFace = false;
doZColorEdge = false;

%% check given handle
try
	findobj(hdl);
catch
	error(['Given handle for drawing XYPlotObject is invalid. ' ...
		'Type ''help redraw'' for more information.']);
end

%% Then plot
axes(hdl);

% give the user the possibility to just specify the startColor; if this is
% what happened, all values should be equal and whole
if isscalar(XYPlotObject.params.ML.Color{1}) && ...
		isnumeric(XYPlotObject.params.ML.Color{1})
	colors = [XYPlotObject.params.ML.Color{:}];
	if ~isequal( ones(1,length(colors)), colors/colors(1) )
		error('Internal error: invalid color specification.')
	end
	startColor = colors(1);
end

CO = get(gca, 'ColorOrder');
XYPlotFieldNames = fieldnames(XYPlotObject.params.ML);

% Run through the plot parameters and build a string for the plot command
fieldArray = cell(1, size(XYPlotFieldNames, 1));
for FNCounter = 1:size(XYPlotFieldNames, 1)
	fieldArray{FNCounter} = XYPlotObject.params.ML.(XYPlotFieldNames{FNCounter});
end

% Set options for line() and do the plotting 
for i = 1:size(XYPlotObject.XData, 1) % each row in XData contains a new plot
	paramStruct = [];
	for FNCounter = 1:size(XYPlotFieldNames, 1)
		paramStruct.(XYPlotFieldNames{FNCounter}) = fieldArray{FNCounter}{i};
	end
	
	%if colors are default, look at the argument startColor, it is used by
	%KPlotObject to give different plotObjects different colors:
	% Set params.ML.Color to the color derived from StartColor
	% params.ML.Color will be further used to set MarkerFaceColor and
	% MarkerEdgeColor if they are still the default values ( = empty).
	% All values: Color, MarkerFaceColor and MarkerEdgeColor are saved here
	% and re-used if the XPobject is redrawn.
	if isequal( '', XYPlotObject.params.ML.Color{i} ) ...
			|| (isscalar( XYPlotObject.params.ML.Color{i} ) && ...
			isnumeric( XYPlotObject.params.ML.Color{i} ))
		XYPlotObject.params.ML.Color{i} = ...
			CO( mod(i + startColor - 1, size(CO, 1)) + 1, :);
		% also change the color in paramStruct
		paramStruct.Color = XYPlotObject.params.ML.Color{i};
	end
	if isequal( 'fill', XYPlotObject.params.ML.MarkerFaceColor{i})
		XYPlotObject.params.ML.MarkerFaceColor{i} = XYPlotObject.params.ML.Color{i};
		paramStruct.MarkerFaceColor = XYPlotObject.params.ML.MarkerFaceColor{i};
		
		%by Abel: set ZColor if needed
		% if the MarkerFaceColor is set to 'fill' we need to include the
		% MarkerFaceColor in Z-Axes coloring if the XYPlotObject.ZColors is
		% not empty
		doZColorFace = true;
	end
	if isequal( '', XYPlotObject.params.ML.MarkerEdgeColor{i})
		XYPlotObject.params.ML.MarkerEdgeColor{i} = XYPlotObject.params.ML.Color{i};
		paramStruct.MarkerEdgeColor = XYPlotObject.params.ML.MarkerEdgeColor{i};
		
		%by Abel: set ZColor if needed
		% See MarkerFaceColor above
		doZColorEdge = true;
	end
	
	%by Abel: implement ZColor
	% - Since matlab can't color single points, each point has to be drawn
	%   seperately using the line() function.
	% - ZColor plotting is invoked automatically when XYPlotObject.ZColors
	%   is set by setZColors() 
	if ~all(strcmpi('none', XYPlotObject.ZColors))
		%Get the data 
		XData = XYPlotObject.XData{i};
		YData = XYPlotObject.YData{i};
		
		%doZColorEdge == doZColorFace == 0
		%This situation may occur if the XYobject was drawn first without
		%ZColoring (which sets face and edge colors to ML.Color) and now 
		%again with ZColors. 
		if ~doZColorEdge && ~doZColorFace
			doZColorEdge = true();
			doZColorFace = ~strcmpi( XYPlotObject.params.ML.MarkerFaceColor{i}, 'none');
		end
		
		%Remove Edge and Face if needed
		% We need to set the MarkerEdgeColor and MarkerFaceColor to their
		% default values again in order to invoke doZColorEdge and
		% doZColorFace when redraw() is executed on the object
		if doZColorEdge
			XYPlotObject.params.ML.MarkerEdgeColor{i} = '';
		end
		if doZColorFace
			XYPlotObject.params.ML.MarkerFaceColor{i} = 'fill';
		end

		%Plot points seperately
		for n=1:length(YData)
			%Color the edge
			if doZColorEdge
				paramStruct.MarkerEdgeColor = XYPlotObject.ZColors{i,n};
			end
			%Color the face
			if doZColorFace
				paramStruct.MarkerFaceColor = XYPlotObject.ZColors{i,n};
			end
			%run line and save handles 
			plotHdl(n) = line(XData(n), YData(n), paramStruct);
		end
		
		%save only the first hdl in .plotHdl for compatibility
 		XYPlotObject.plotHdl(i) = plotHdl(1);
		%save all in other place
		XYPlotObject.pointHdl{i} = plotHdl;
	else
		XYPlotObject.plotHdl(i) = ...
			line(XYPlotObject.XData{i}, XYPlotObject.YData{i}, paramStruct);
	end
end
