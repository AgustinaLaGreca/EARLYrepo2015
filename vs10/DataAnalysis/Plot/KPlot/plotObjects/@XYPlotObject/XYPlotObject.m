function newXYPlot = XYPlotObject(XData, YData, varargin)
% XYPLOTOBJECT Creates an XYPlotObject instance.
%
% An XYPlotObject contains data and information for a rasterplot or a
% vectorplot. Plotting happens using the redraw method. After construction
% the XYPlotObject can be added to a Panel object.
%
% newXYPlot = XYPlotObject(XData, YData, paramList) 
% Creates an XYPlotObject object from the given vectors. Parameter pairs
% can be given in paramList.
%
%       XData: XData for the plots. It is a column cell array containing a
%              list of row vectors, or a two dimensional numeric array. In
%              the latter case, all vectors need to have the same length.
%       YData: YData for the plots. It is a column cell array containing a
%              list of row vectors, or a two dimensional numeric array. In
%              the latter case, all vectors need to have the same length.
%      Marker: Markers for the plots; can be a single value (same for every
%              plot) or a column cell array  with a value for every plot.
%       Color: Colors for the plots; can be a single value (same for every
%              plot) or a column cell array  with a value for every plot.
%MarkerFaceColor: An extra mode was added for MarkerFaceColor, the mode 'fill'
%                 sets the MarkerFaceColor to the value of Color
%   LineStyle: Linestyle for the plots; can be a single value (same for 
%              every plot) or a column cell array  with a value for every
%              plot.
% 
% When called without any arguments, an empty XYPlotObject is returned,
% this can be used to preallocate arrays.
% 
% EXAMPLE
%	%Draw a line and a single dot
%	X{1,:} = [1 2 3];	%X-Values of line
%	Y{1,:} = [1 2 3];	%Y-Values of line
%	X{2,:} = [4];		%X-Value of dot
%	Y{2,:} = [4];		%Y-Value of dot
%	Markers{1,:} =	['none'];	%No markers needed for the line
%	Markers{2,:} =	['o'];		%Use 'o' as marker for the dot
%	Obj = XYPlotObject(X, Y, 'marker', Markers)

% Created by: Kevin Spiritus
% Edited by: Ramses de Norre

%% ---------------- CHANGELOG -----------------------
%  Fri Jan 28 2011  Abel   
%	- Added example in Doc
%  Tue Apr 26 2011  Abel   
%   - Reorganised code so XYPlotObject() without arguments returns an
%   object with default parameters.
%  Fri Aug 12 2011  Abel   
%   - Added options for Z-Axes coloring
%  Thu Aug 18 2011  Abel   
%   - Added ZNaNColor option


%% ---------------- TODO ----------------------------
%  Thu Oct 27 2011  Abel   
%  1) Add ZData to the object:
%	- Follow the same flowchart as for X/Y data 
%   -> When Z is empty we use YData
% 
%% ---------------- CHANGELOG -----------------------
%  Tue Nov 22 2011  Abel   
%   - Implemented TODO nr 1

%% standard properties
% Internals
newXYPlot.XData						= [];
newXYPlot.YData						= [];
newXYPlot.plotHdl					= [];
%by Abel: save Hdls of single points if they were drawn separately (for ex. by
%color-by-Z-axes plotting)
newXYPlot.pointHdl					= [];

% User adaptable MatLab options
newXYPlot.params.ML.Marker			= {'x'};
newXYPlot.params.ML.MarkerEdgeColor = {''};
newXYPlot.params.ML.MarkerFaceColor = {'none'};
newXYPlot.params.ML.MarkerSize		= {6};
newXYPlot.params.ML.LineStyle		= {'-'};
newXYPlot.params.ML.LineWidth		= {0.5};
newXYPlot.params.ML.Color			= {''};
newXYPlot.params.ML.ButtonDownFcn	= {''};
newXYPlot.params.ML.Visible         = {'on'};

% by Abel: new color-by-Z-axes option
newXYPlot.ZColors = {'none'};
newXYPlot.ZScale = [];
newXYPlot.ZData = [];
newXYPlot.params.ZNaNColor = {'none'};

% Additional parameters (KPlot, not matlab)
newXYPlot.params.RedrawOnResize = false;

% internal vars 
gotZ = false();
ZData = [];

%% Parameter checking
% Return empty objects with default values if no arguments (data) was given
if isequal(nargin,0)
    newXYPlot = class(newXYPlot, 'XYPlotObject');
    return
end
% Check min Nr of Data 
if nargin < 2
    error(['XYPlotObject expects at least two arguments or none. ' ...
        'Type ''help XYPlotObject'' for more information.']);
end
% Test if ZData was given 
if nargin >= 3
	if iscell(varargin{1})
		gotZ = all(cellfun(@isnumeric, varargin{1}));
		if all(cellfun(@isempty, varargin{1}))
			varargin{1} = [];
		end
	else
		gotZ = isnumeric(varargin{1});
	end
end
% Save parameter list and Zdata if needed 
if gotZ
	ZData = varargin{1};
	varargin = varargin(2:end);
else
	ZData = YData;
end
% If ZData is empty, use YData
if isempty(ZData)
	ZData = YData;
end
% Save rest of parameters, these will be checked later, by processParams
paramsIn = varargin; 


% If XData is simply a list of numeric values, check X/Y consistency and
% save their values in a cell 
if isnumeric(XData)
	%Check X/Y consistency
	if ~isequal( 2, ndims(XData) ) || ~isequal(size(XData), size(YData), size(ZData))  || ...
            ~isnumeric(YData) || ~isnumeric(ZData)
        error(['Wrong format of the arguments. ' ...
            'Type ''help XYPlotObject'' for more information.']);
	end
	
    %Each row is a new line on the plot => save X/Y in cell with X/Y values 
	%saved per row.
	length = size(XData, 1);
    XDataTemp = cell(1, length);
    YDataTemp = cell(1, length);
	ZDataTemp = cell(1, length);
    for i = 1 : length
        XDataTemp{i} = XData(i, :);
        YDataTemp{i} = YData(i, :);
		ZDataTemp{i} = ZData(i, :);
    end
    XData = XDataTemp';
    YData = YDataTemp';
	ZData = ZDataTemp';
end

% If XData and YData was not a list of values, they should be cells
if ~iscell(XData) || ~iscell(YData) || ~iscell(ZData)
    error(['Wrong format of the arguments. ' ...
        'Type ''help XYPlotObject'' for more information.']);
end

% Check X/Y length consistency
if ~isequal( 2, ndims(XData) ) || ~isequal( 1, size(XData,2) ) || ...
        ~isequal(size(XData), size(YData), size(ZData))
    error(['Wrong format of the arguments. ' ...
        'Type ''help XYPlotObject'' for more information.']);
end

% Delete empty rows
i = 1;
while i <= size(XData,1)
    if isempty(XData{i})
        XData = {XData{[(1:i-1), (i+1:end)]}}';
        YData = {YData{[(1:i-1), (i+1:end)]}}';
		ZData = {ZData{[(1:i-1), (i+1:end)]}}';
    else
        i = i+1;
    end
end

% Continue checking data format: XData & YData should now be cells of equal size 
% with numeric values and a row per plotted line
for i=1:size(XData,1)
    if ~isnumeric(XData{i}) || ~isnumeric(YData{i}) || ~isnumeric(ZData{i})
        error(['Wrong format of the arguments. ' ...
            'Type ''help XYPlotObject'' for more information.']);
    end
    if ~isequal( 2, ndims(XData{i}) ) || ~isequal( 1, size(XData{i},1) ) || ...
            ~isequal(size(XData{i}), size(YData{i}), size(ZData{i}))
        error(['Wrong format of the arguments. ' ...
            'Type ''help XYPlotObject'' for more information.']);
    end
end

%% Assign data
newXYPlot.XData = XData;
newXYPlot.YData = YData;
newXYPlot.ZData = ZData;

%% Handle parameters
nRows = size(XData, 1); % The amount of vectors in input data

% define structure of parameters: just fill with defaults based on the
% number of input XData
optFields = fieldnames(newXYPlot.params.ML);
for n=1:size(optFields, 1)
	newXYPlot.params.ML.(optFields{n}) = repmat(newXYPlot.params.ML.(optFields{n}), nRows, 1);
end

% Parameters given as arguments are placed in the params structure.
% Other entries remain 'def'. This def value is replaced when plotting.
newXYPlot.params = processParams(paramsIn, newXYPlot.params);

%% Return object
newXYPlot = class(newXYPlot, 'XYPlotObject');
