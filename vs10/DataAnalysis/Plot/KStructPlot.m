function [argOut, params] = KStructPlot(varargin)
% KStructPlot - Create a plot using the fields of a structure-array
%   KStructPlot(S1, XFieldNames1, YFieldNames1, S2, XFieldNames2, YFieldNames2,
%               ..., parameters)
%
%   The data for the plots is contained in the structures S1, S2, .... The
%   fields containing X and Y data are indicated in XFieldNames and
%   YFieldNames parameters.
%
%   Each row of these data columns can contain a scalar, represented by a
%   dot in the plot, or a vector, represented by a line.
%
%   Different structures get different colors in the plot.
%
%   Parameters:
%         colors          : the colors used for different structures
%         linestyles      : the linestyles used for different structures
%         markers         : the markers used for different structures,
%                           appending an 'f' to a marker symbol makes it filled
%		  markersizes	  : MarkerSize for points
%         execevalfnc     : Matlab statement to execute when clicking on a
%                           datapoint or curve. Fieldnames in these statements
%                           must be enclosed between dollar signs and for
%                           branched structures fieldnames can be given using
%                           the dot as a fieldname separator
%         animalidfield   : the name of the field that identifies the animal for
%                           a given structure-array row. A cell-array of strings
%                           can be supplied if this field is different for
%                           different structure-arrays
%         cellidfield     : the name of the field containing cell id's
%         totalidfields   : TODO
%         dispstats       : display statistics (animals, cells, ...)
%         indexexpr       : condition for including an index
%         intersectfields : only include indexes from the intersection of these
%                           fields
%         info            : TODO
%         xlim            : two element vector giving the lower and upper x
%                           limits
%         ylim            : two element vector giving the lower and upper y
%                           limits
%         gutter          : boolean, whether or not to display a gutter
%                           with points that are not inside the plot range
%         johnson         : yes/[no], whether or not to plot in the Johnson
%                           scale
%     ignoreTotalIDFields : yes/[no], whether we should try to deduce the
%                           totalidfields or not
%		  returnstruct	  : return plotted data in struct
%		  fitfnc          :	fit a cubic spline function through data
%							'none' or 'spline'
%		  fitsampleratio  : Ratio which defines oversampling of fit function
%							to be plotted.
%		  title			  : Plot title
%		  xlabel/ylabel   : Axis labels
%         logX/logY       : Axis in log scale (yes/no)
%         logXformat/logYformat :	Useful in situations where the X-as is in log scale.
%									Set to 'auto' for matlab labels or leave blank to use
%									xlog125 as labels
%		  xTickLabels     : Print a user defined string on the xTicks
%         legend          : Add a user defined legend instead of the struct
%                           names. You can add a value from a column of the
%                           input struct by using the dollar sign syntax.
%         legendParam     : Change the default options for the legend
%
%         showexeceval    : (false/true) Show ButtonDownFcn for debugging
%   Run Kstructplot() without arguments to see default values for these
%   parameters.
%	OUTPUT:
%		argOut:	Struct containing plotted values in the same format as
%					the input struct array. If returnstruct == true
%
%
% SEE ALSO panel legend


%% ---------------- CHANGELOG -----------------------
%[Thu Jan 15 2015 (Abel)]: Bugfix: names of input variables are now more
%reasonably handled. This renders support for D(1) type of inputs 
%
%  Mon Jan 31 2011  Abel
%   bugfix in localGetButtonDownFcn: always row vector: '' -> should be .'
%	bugfix in localGetButtonDownFcn: param.totalidfields is cell array
%  Tue Feb 1 2011  Abel
%   bugfix: ignoreTotalIDFields = 'yes' was not checked in localGetButtonDownFcn
%  Thu Feb 3 2011  Abel
%   Added function output:
%		argOut: Struct array containing plotted values in the same format as
%		the input struct array. If returnstruct == true
%  Fri Feb 4 2011  Abel
%   bugfix in output struct
%	added support for by-inputstruct-execevalfnc
%  Tue Feb 8 2011  Abel
%   - Added test for number of input elements in a cell input parameter
%  Thu Feb 10 2011  Abel
%   - bugfix output struct when no input name can be retrieved
%  Thu Feb 10 2011  Abel
%   - Added markersize option
%	- Added spline fitting
%  Mon Feb 28 2011  Abel
%   - Added extra plotting options: title, x&y-label,..etc.
%  Fri Mar 4 2011  Abel
%   - Added option to keep standard matlab log-scale
%  Fri Apr 8 2011  Abel
%   - Added reverseX/Y option
%  Thu Apr 14 2011  Abel
%   - Bugfix in totalidfield filename retrieval
%  Tue Apr 19 2011  Abel
%   - Added legendObject support
%  Fri Aug 12 2011  Abel
%   - Added ZColor support
%   - Options in defaultParams.checkfornoe are now repmat()'ed to the correct
%     number of values when only one option was given.
%   - Added data.idx so we can see the selections that have been done on
%     data.Y
%  Thu Aug 18 2011  Abel
%   - Added support for ZNaNColor option
%  Mon Aug 22 2011  Abel
%   - Added showexeceval option
%  Wed Sep 7 2011  Abel
%   - Changed default legendParams
%  Tue Nov 22 2011  Abel   
%   - Partial rewrite: Actual plotting was shifted to KStructPlotPage()
%   - Partial rewrite: Gutters are within separate Panel() Objects.
%  Tue Apr 24 2012  Abel   
%   - added support for edge histograms
%   - bugfix: evalfnc did take filtered out X/Y values into account
%[Thu Jan 29 2015 (Abel)]: 
%- Added support for dollar replacement syntaxt in legend function
%[Mon Mar 9 2015 (Abel)]: 
% - Added support for parameters in sturct for (paramstruct option)
% - Added totalidfields to the options to check for the number of inputs 
% - Allow the user to override an incorrect nr op options 
%[Thu Mar 19 2015 (Abel)]:
% - bugfix totalidfields: when fieldname is not "filename"...


%% ---------------- TODO -----------------------
%  Fri Aug 12 2011  Abel
%  1) Check for each defaultParams.checkfornoe if their nr of elements is
%  re-checked needlessly locally to their code-function
%[Thu Jan 29 2015 (Abel)]: 
% - check gutter, always circles on edges? 
% - check possible errors when parsing totalID field 

%% Defaults
defaultParams.colors              = {'b'; 'g'; 'r'; 'c'; 'm'; 'y'};
defaultParams.linestyles          = {'-'};
defaultParams.markers             = ...
	{'o'; '^'; 'v'; '*'; '+'; 'x'; '<'; '>'; 's'; 'd'; 'p'; 'h'};
defaultParams.markersizes		  = {6};
defaultParams.execevalfnc         = '';
defaultParams.animalidfield       = {};
defaultParams.cellidfield         = {};
defaultParams.totalidfields       = {'ds.filename', 'ds.seqid'};
defaultParams.dispstats           = 'yes';
defaultParams.indexExpr           = '';
defaultParams.intersectfields     = {'',''};
defaultParams.info                = 'yes';
defaultParams.xlim                = [-inf, inf];
defaultParams.ylim                = [-inf, inf];
defaultParams.gutter              = false;
defaultParams.johnson             = 'no';
defaultParams.ignoreTotalIDFields = 'no';
defaultParams.returnstruct        = false;
defaultParams.fitfnc			  = {'none'};
defaultParams.fitsampleratio	  = 5;
defaultParams.title		          = '';
defaultParams.xlabel			  = '';
defaultParams.ylabel			  = '';
defaultParams.logX				  = '';
defaultParams.logY				  = '';
defaultParams.xTickLabels		  = '';
defaultParams.logXformat		  = '';
defaultParams.logYformat		  = '';
defaultParams.reverseX			  = '';
defaultParams.reverseY            = '';
defaultParams.legend              = '';
defaultParams.legendParam         = {'Location', 'NorthEast'};
defaultParams.ticksdir            = 'out';

% by abel: Parameters in struct support
defaultParams.paramstruct = '';       %struct to read parameters from. These options will be 
                                      %overwritten by any options you specify directly on the 
                                      %commandline
% by abel: when the number of options isn't correct, do we produce an error
% or a warning? 
defaultParams.strictParamCount = true;

% Z-Axes color plotting
defaultParams.ZColor = 'none';
defaultParams.ZColorMap = 'jet';
defaultParams.ZScale = [];
defaultParams.ZData = '';
defaultParams.ZColorType = 'fill';
defaultParams.ZNaNColor = 'none';
defaultParams.addZColorbar = false;

% Edge histogram
defaultParams.hist = false;
defaultParams.histbinnr = 128;
defaultParams.histyunit = '#';


% Debug execeval function
defaultParams.showexeceval = false; %Show ButtonDownFcn for debugging

% Variables used in script
defaultParams.checkfornoe = {'execevalfnc', 'indexExpr', 'legend', 'ZData',...
	'ZScale', 'totalidfields'};
defaultParams.plotoptions = {'title',...
	'xlabel',...
	'ylabel',...
	'logX',...
	'logY',...
	'xTickLabels',...
	'logXformat',...
	'logYformat',...
	'reverseX',...
	'reverseY',...
	'addZColorbar',...
	'ticksdir'};            %Panel() property (TickDir in Axes props)

% Return factory defaults?
argOut = []; % standard, return nothing
if nargin < 1
	argOut = defaultParams;
	return
end

%% main function
% parse input
%by Abel:
% Save the KStructplot input as a string in order to show it on the plot
% lateron.
% - get the names of the input variables
% - build commandline input as string in inputString
% The inputnames are saved in inputNames which is used in localLoadData()
% to prepare a legend on the plot. length(inputNames) must be ==
% length(varargin), although the names for string inputs don't matter since
% we only look into inputNames for arguments we know are structs.
inputNames = cell(1, nargin);
inputString = 'KStructplot(';
for cArgin = 1:nargin
    %by Abel:
    % Attention! an input struct array D, can only be resolved by inputname
    % if it is given on the commandline as D and not as one element D(1) for example. 
    % String inputs like 'curve.indepval' are also not resolvable by
    % inputname()
	inputNames{cArgin} = inputname(cArgin);
	if isempty(inputNames{cArgin})
        if isa(varargin{cArgin}, 'struct')
            inputNames{cArgin} = sprintf('Unknown_%d', cArgin);
        elseif isa(varargin{cArgin}, 'cell')
            inputNames{cArgin} = cell2str(varargin(cArgin));
        elseif isa(varargin{cArgin}, 'char')
            inputNames{cArgin} = sprintf('''%s''', varargin{cArgin});
        elseif isa(varargin{cArgin}, 'logical')
            inputNames{cArgin} = sprintf('%d', varargin{cArgin});
        elseif isnumeric(varargin{cArgin})
            inputNames{cArgin} = sprintf('%d', varargin{cArgin});
        else
            inputNames{cArgin} = sprintf('Unknown_input_type_%d', cArgin);
        end
        
        % add to varString
        varString = inputNames{cArgin};
		
	else
		varString = inputNames{cArgin};
	end
	
	if cArgin > 1 && cArgin < nargin
		inputString = strcat(inputString, ',' , varString);
	else
		inputString = [inputString varString];
	end
end
inputString = [inputString, ');'];

% parse input params and load data
[data, params, stats] = localLoadData(varargin, inputNames, defaultParams);
dataLength = length(data);

% initialize plot objects
XYPlotLines = repmat(XYPlotObject(), 1, dataLength);
XYPlotPoints = repmat(XYPlotObject(), 1, dataLength);

% declare some leafs
params.xmax = [];
params.xmin = [];
params.ymax = [];
params.ymin = [];
params.alldata = struct('x', [], 'y', [], 'z', []);


% start loop over input structures
for cStruct = 1:dataLength
	% collect aditional info
	sortedCellData = localSortCellData(data(cStruct), params);
	
	stats(cStruct).nDots = 0;
	for cRow = 1:length(data(cStruct).rows)
		X = data(cStruct).rows(cRow).X;
		Y = data(cStruct).rows(cRow).Y;
		% Count only the points that are inside the given x and y limits
		stats(cStruct).nDots = stats(cStruct).nDots + ...
			length(X( X >= params.xlim(1) & X <= params.xlim(2) & ...
			Y >= params.ylim(1) & Y <= params.ylim(2) ));
	end
	
	stats(cStruct).nCells = 0;
	for cRow = 1:length(sortedCellData)
		X = sortedCellData(cRow).X;
		Y = sortedCellData(cRow).Y;
		% Count only the cells that have points inside the given x and y limits
		if length(X( X >= params.xlim(1) & X <= params.xlim(2) & ...
				Y >= params.ylim(1) & Y <= params.ylim(2) )) >= 1
			stats(cStruct).nCells = stats(cStruct).nCells + 1;
		end
	end
	
	if params.gutter
		params.gutterLimits = calcGutterLimits(params.xlim, params.ylim);
	else
		params.gutterLimits = {params.xlim, params.ylim};
	end
	
	% Create plot objects
	colorIdx = rotIdx( cStruct, length(params.colors) );
	markerIdx = rotIdx( cStruct, length(params.markers) );
	lineIdx = rotIdx( cStruct, length(params.linestyles) );
	markersizeIdx = rotIdx( cStruct, length(params.markersizes) );
	fitfncIdx = rotIdx( cStruct, length(params.fitfnc) );
	
	% We create 2 layers:
	% - one layer for the connecting lines
	% - one layer for the dots (with point & click functionality)
	
	sortedX = {sortedCellData.X}';
	sortedY = {sortedCellData.Y}';
	
	% Adding an 'f' to a Marker symbol gives it the color of the border
	if length(params.markers{markerIdx}) > 1 && ...
			isequal(params.markers{markerIdx}(2), 'f')
		markerFaceColor = params.colors{colorIdx};
	else
		markerFaceColor = 'none';
	end
	
	X = {data(cStruct).rows(:).X}';
	Y = {data(cStruct).rows(:).Y}';
    
	%by Abel: Z-Axes plotting for dots
	Z = {data(cStruct).rows(:).Z}'; 
	if all(cellfun(@isempty, Z))
		Z = Y;
	end
	
	% Inverting of Y data is done here, the other adjustments for the
	% Johnson scale are done in localPlotPanel
	if isequal(params.johnson, 'yes')
		Y = cellfun(@(y) 1 - y, Y, 'UniformOutput', false);
		sortedY = cellfun(@(y) 1 - y, sortedY, 'UniformOutput', false);
	end
	
	if params.gutter
		[X, Y, Z, gutter(cStruct)] = createGutter(X, Y, Z, params.xlim, params.ylim, ...
			params.gutterLimits);
	end
	
	%Remove any empty Y entries
	notEmptyIdx = ~cellfun(@isempty, Y);
	X = X( notEmptyIdx );
	Y = Y( notEmptyIdx );
	Z = Z( notEmptyIdx );
	
	%by Abel: fit cubic spline
	if strcmp(params.fitfnc{fitfncIdx}, 'spline')
		[sortedX, sortedY] = fitSpline_(sortedX, sortedY, params.fitsampleratio);
	end
	
	% dots
	XYPlotPoints(cStruct) = XYPlotObject(X, Y, Z,...
		'color', params.colors{colorIdx}, ...
		'marker', params.markers{markerIdx}(1), ...
		'markerfacecolor', markerFaceColor, ...
		'linestyle', 'none', ...
		'ButtonDownFcn', {data(cStruct).rows(notEmptyIdx).buttonDownFcn}.',...
		'MarkerSize', params.markersizes{markersizeIdx} );
	
	% Connecting lines
	XYPlotLines(cStruct) = XYPlotObject(sortedX, sortedY, ...
		'color', params.colors{ colorIdx }, 'marker', 'none', ...
		'linestyle', params.linestyles{ lineIdx });
    
	% By Abel: return the plotted data in the original format of the input
	% struct. Gutter should not be included
	if (params.returnstruct)
		name = data(cStruct).name;
		if isfield(argOut, name)
			name = sprintf('%s_%d', data(cStruct).name, cStruct);
		end
		argOut.(name) = returnPlottedData(X, Y, data(cStruct), params.xlim, params.ylim);
	end
	
	% by Abel: save max/min values for later
	params.xmax = max([X{:} params.xmax]);
	params.xmin = min([X{:} params.xmin]);
	params.ymax = max([Y{:} params.ymax]);
	params.ymin = min([Y{:} params.ymin]);
	
	% save combine X/Y/Z data for later
 	params.alldata.x = [params.alldata.x, X{:}];
 	params.alldata.y = [params.alldata.y, Y{:}];
 	params.alldata.z = [params.alldata.z, Z{:}];
	
end

%% Create KPlot Panel()s and plot using KStructPlotPage()
pPageParams = struct();
pPageParams.inputstring = inputString;

% - input string for use in datePanel()
params.inputstring = inputString;

% - create panel for main plot
mainPanel = localPlotPanel(XYPlotLines, XYPlotPoints, {data(:).XFieldNames}, ...
	{data(:).YFieldNames}, params);

% - Add statistics statistics
if isequal('y', lower(params.info(1)))
	statsText = localPlotStats({data(:).name}', stats);
	pPageParams.statsstring = statsText;
end

%create gutter panels
params.gutterparam = [];
pPageParams.gutterpaneltop = [];
pPageParams.gutterpanelbottom = [];
pPageParams.gutterpanelleft = [];
pPageParams.gutterpanelright = [];
if params.gutter
	% Save gutter
	params.gutterparam = gutter;
	[ pPageParams.gutterpaneltop, pPageParams.gutterpanelbottom,...
		pPageParams.gutterpanelleft, pPageParams.gutterpanelright ]  = createGutterPanels( params );
end

%create edge histograms
pPageParams.ehistx = [];
pPageParams.ehisty = [];
if params.hist
	[pPageParams.ehistx, pPageParams.ehisty] = createHistPanels( params ); 
end

% - create plot using KStrunctPlotPage()
pPageParams.mainpanel = mainPanel;
%save struct names for writing legends
params.structnames = {data(:).name};
[ axH, param ] = KStructPlotPage(params, pPageParams);


%% createGutter
% Function executed seperately for each input struct
function [X, Y, Z, gutter] = createGutter(X, Y, Z, xlim, ylim, gutterLimits)

[lGutX, lGutY, rGutX, rGutY, tGutX, tGutY, bGutX, bGutY] = deal({});


gutter.top = [];
gutter.bottom = [];
gutter.left = [];
gutter.right = [];

xFiltData = {};
yFiltData = {};
zFiltData = {};

for n = 1:length(X)
	leftXidx = ([X{n,:}] <= xlim(1));
	rightXidx = ([X{n,:}] >= xlim(2));
	bottomYidx = ([Y{n,:}] <= ylim(1));
	topYidx = ([Y{n,:}] >= ylim(2));
	
	xData = X{n,:};
	yData = Y{n,:};
	zData = Z{n,:};
	
	if any(leftXidx ~= 0)
		%collect gutter
		lGutX{end+1} = zeros(1, length(leftXidx)); %empty rows will be removed by XYPlotObject()
		lGutX{end} = lGutX{end}(leftXidx);
		lGutY{end+1} = Y{n,:}(leftXidx);
		%filter X/Y data
		xData = xData( not(leftXidx) );
		yData = yData( not(leftXidx) );
		zData = zData( not(leftXidx) );
	elseif any(rightXidx ~= 0)
		%collect gutter
		rGutX{end+1} = zeros(1, length(rightXidx)); %empty rows will be removed by XYPlotObject()
		rGutX{end} = rGutX{end}(rightXidx);
		rGutY{end+1} = Y{n,:}(rightXidx);
		%filter X/Y data
		xData = xData( not(rightXidx) );
		yData = yData( not(rightXidx) );
		zData = zData( not(rightXidx) );
	elseif any(bottomYidx ~= 0)
		%collect gutter
		bGutY{end+1} = zeros(1, length(bottomYidx));
		bGutY{end} = bGutY{end}(bottomYidx);
		bGutX{end+1} = X{n,:}(bottomYidx);
		%filter X/Y data
		xData = xData( not(bottomYidx) );
		yData = yData( not(bottomYidx) );
		zData = zData( not(bottomYidx) );
	elseif any(topYidx ~= 0)
		%collect gutter
		tGutY{end+1} = zeros(1, length(topYidx));
		tGutY{end} = tGutY{end}(topYidx);
		tGutX{end+1} = X{n,:}(topYidx);
		%filter X/Y data
		xData = xData( not(topYidx) );
		yData = yData( not(topYidx) );
		zData = zData( not(topYidx) );
	end
	
	%save if not empty
	xFiltData{end+1} = xData;
	yFiltData{end+1} = yData;
	zFiltData{end+1} = zData;

end

%Filtered Data 
X = tocol(xFiltData);
Y = tocol(yFiltData);
Z = tocol(zFiltData);

%Gutter
gutter.top = { tGutX, tGutY };
gutter.bottom = { bGutX, bGutY };
gutter.left = { lGutX, lGutY };
gutter.right = { rGutX, rGutY };
%% gutterLimits
function gutterLimits = calcGutterLimits(xlim, ylim)
xglim = xlim;
yglim = ylim;
xdiff = abs(diff(xlim))*.02;
ydiff = abs(diff(ylim))*.02;

xglim(1) = xglim(1) - xdiff;
yglim(1) = yglim(1) - ydiff;
xglim(2) = xglim(2) + xdiff;
yglim(2) = yglim(2) + ydiff;

gutterLimits = { xglim, yglim};
%%

%% localLoadData
function [data, params, stats] = localLoadData(varargin, inputNames, defaultParams)
% Get Structures, fieldnames and extra arguments from the command line
% parameters.

% Parse varargin
structIn = {};    % the structures themselves
XFieldNames = {}; % fieldnames for X data
YFieldNames = {}; % fieldnames for Y data
cStruct = 0;      % number of handled plot
while ~isempty(varargin) && isstruct(varargin{1})
	cStruct = cStruct + 1;
	% Format is: KStructPlot(S1, XFiendNames1, YFieldNames1, S2,
	% XFiendNames2, YFieldNames2, ..., parameters)
	if ~ischar(varargin{2}) || ~ischar(varargin{3})
		error(['Wrong arguments for ' mfilename '.']);
	end
	structIn{cStruct} = varargin{1};
	XFieldNames{cStruct} = varargin{2};
	YFieldNames{cStruct} = varargin{3};
    
	% move on
	%by abel: replace { {} } by ()
	%varargin = {varargin{4:end}};
	varargin = varargin(4:end);
end

stringsIdx = find(cellfun(@(s) isa(s, 'char'), varargin));
strings = varargin(stringsIdx);
% Determine filename if it is not set manually, we look for options ds, ds1
% and ds2
if ~ismember('totalidfields', strings)
	defaultParams.totalidfields = {};
	for cStruct = 1:length(structIn)
		S = structIn{cStruct};
		possibilities = {{'ds.filename', 'ds.seqid'}, ...
			{'ds1.filename', 'ds1.seqid'}, {'ds2.filename', 'ds2.seqid'}, ...
            {'Ds.filename', 'Ds.seqid'}};
		possibilityIdx = [possibilities{isfield(S, {'ds', 'ds1', 'ds2','Ds'})}];
		if ~isempty(possibilityIdx)
			defaultParams.totalidfields{cStruct} = possibilityIdx;
		end
	end
end

%by abel:
% add support for input parameters in struct format
[ gotParamStruct, paramIdx ] = ismember('paramstruct', strings);
if gotParamStruct
    paramIdx = stringsIdx(paramIdx);
    defaultParams = processParams( varargin{paramIdx+1}, defaultParams);
    if length(varargin) == paramIdx+1
        varargin = varargin(1:paramIdx-1);
    else
        varargin = varargin([1:paramIdx-1, paramIdx+2:end]);
    end
end
params = processParams(varargin, defaultParams);

%By Abel: Check if the number of inputs in cells equals the number of
%structs
% - a single input is automatically duplicated to the number of structs
for noe = 1:length(params.checkfornoe)
	isCell = iscell(params.(params.checkfornoe{noe}));
	isEmpty = isempty(params.(params.checkfornoe{noe}));
	
	%If we get a single input which is not a cell -> put it in a cell
	if ~isCell && ~isEmpty
		params.(params.checkfornoe{noe}) = {params.(params.checkfornoe{noe})};
	end
	
	nCell = length(params.(params.checkfornoe{noe}));
	nStruct = length(structIn);

	%If there is only one element -> repmat to the number of input structs
	if nCell == 1
		params.(params.checkfornoe{noe}) = repmat(params.(params.checkfornoe{noe}),1,nStruct);
		%If by now, the nr ~= nr of structs. Die with an error
	elseif isCell && nCell ~= nStruct
        
        doDie = params.strictParamCount;
        if doDie
            error('The number of cell elements in option:%s is not equal to the number of input structs:%d', params.checkfornoe{noe}, length(structIn));
        else
            warning('The number of cell elements in option:%s is not equal to the number of input structs:%d', params.checkfornoe{noe}, length(structIn));
        end

	end
end

%by Abel: Add legend
% Add legend
% - as input for the legend option: { 'a cell string of names, length = nr
% of structs' }
% - a colum of the input struct {'$fieldname$)
% if a specific legend was supplied, check to see if we need to resolve
% a fieldname
if ~isempty(params.legend)
    for cStruct = 1:length(structIn)
        if any(strfind(params.legend{cStruct}, '$'))
            myInStruct = structIn{cStruct};
            legendStr = eval( parseExpression(params.legend{cStruct}, 'myInStruct(1)'));
            params.legend{cStruct} = any2str(legendStr);
        end
    end
end
    


% load the data for the structs
data = [];
stats = [];
for cStruct = 1:length(structIn)
	data(cStruct).rows = localGetStructData(structIn{cStruct}, ...
		XFieldNames{cStruct}, YFieldNames{cStruct}, params, cStruct);
	data(cStruct).XFieldNames = XFieldNames{cStruct};
	data(cStruct).YFieldNames = YFieldNames{cStruct};
	data(cStruct).name = inputNames{(cStruct-1)*3+1};
end

for cStruct = 1:length(structIn)
	if ~isempty(params.animalidfield)
		if ~isequal(length(params.animalidfield), length(structIn))
			error('You must specify as many animalidfields as there are structs.');
		end
		[values, names] = destruct(structIn{cStruct});
		identifiers = {values{:, ...
			cellfun(@(c) isequal(c,params.animalidfield{cStruct}), names)}}; 
		% Only if the numbers in the identifier differ, do they
		% represent different animals
		stats(cStruct).nAnimals = ...
			length(unique(regexprep(unique(identifiers), '[^0-9]', '')));
	else
		stats(cStruct).nAnimals = length(unique(regexprep( ...
			unique({data(cStruct).rows.fileName}), '[^0-9]', '')));
	end
end


%%

%% localGetStructData
function data = localGetStructData(structIn, XFieldNames, YFieldNames, ...
	params, cStruct)
% fill a structure with essential information from structIn
dollarPosX = strfind(XFieldNames, '$');
isDollarXEmpty = isempty(dollarPosX);
dollarPosY = strfind(YFieldNames, '$');
isDollarYEmpty = isempty(dollarPosY);

data.X = [];
data.Y = [];
data.Z = [];
data.idx.X = []; % by Abel, save selected idx for later use
data.idx.Y = [];
data = repmat(data, 1, length(structIn));

%Loop over the input structs
for cRow = 1:length( structIn )
	if isDollarXEmpty
		data(cRow).X = retrieveField( structIn(cRow), XFieldNames );
		data(cRow).X = data(cRow).X{:};
	else
		data(cRow).X = eval(parseExpression( XFieldNames, 'structIn(cRow)' ) );
	end
	
	if isDollarYEmpty
		data(cRow).Y = retrieveField( structIn(cRow), YFieldNames );
		data(cRow).Y = data(cRow).Y{:};
	else
		data(cRow).Y = eval( parseExpression( YFieldNames, 'structIn(cRow)' ) );
	end
	
	% before checking vector lenghts, process intersectfields
	if ~isempty(params.intersectfields{1}) && ~isempty(params.intersectfields{2})
		original{1} = retrieveField(structIn(cRow), params.intersectfields{1});
		original{2} = retrieveField(structIn(cRow), params.intersectfields{2});
		
        original{1} = original{1}{1};
        original{2} = original{2}{1};
        
		if ~isequal(size(original{1}), size(data(cRow).X)) || ...
				~isequal(size(original{2}), size(data(cRow).Y))
			error('Intersect fields should have the same sizes as the data fields.');
		end
		
		I = intersect(original{1}, original{2});
		idx1 = zeros(1, length(I));
		idx2 = zeros(1, length(I));
		for i=1:length(I)
			idx1(i) = find(original{1} == I(i));
			idx2(i) = find(original{2} == I(i));
		end
		% Now limit the data to these indices
		data(cRow).X = data(cRow).X(idx1);
		data(cRow).Y = data(cRow).Y(idx2);
		
		%by Abel: Save selection Idx for later use
		data(cRow).idx.X{end+1} = idx1;
		data(cRow).idx.Y{end+1} = idx2;
	end
	
	% make row vectors
	if ~isequal( 1, size(data(cRow).X, 1) )
		if ~isequal( 1, size(data(cRow).X, 2) )
			error('XData and YData should be vectors or scalars.');
		end
		
		data(cRow).X = data(cRow).X';
		data(cRow).Y = data(cRow).Y';
	end
	
	if ~isequal( size(data(cRow).X), size(data(cRow).Y) )
		error('XData and YData should have the same size');
	end
	
	if ~isempty(params.indexExpr)
		if isa(params.indexExpr, 'cell')
			indexExprOK = find(eval(parseExpression(params.indexExpr{cStruct}, 'structIn(cRow)')));
		else
			indexExprOK = find(eval(parseExpression(params.indexExpr, 'structIn(cRow)')));
		end
		data(cRow).X = data(cRow).X(indexExprOK);
		data(cRow).Y = data(cRow).Y(indexExprOK);
		
		%by Abel: Save selection Idx for later
		data(cRow).idx.X{end+1} = indexExprOK;
		data(cRow).idx.Y{end+1} = indexExprOK;
	end
	
	% determine filename
	%by Abel:
	%Format for totalidfields: structIn = single struct array with identical totalidfields
	% if totalidfields is defined for each input struct array ->
	% totalidfields = param.totalidfields{cStruct} and
	% length(param.totalidfields)  = cStruct since param.totalidfields = {
	% {id1, id2}, ...cStruct }
	% if totalidfields is identical for each input struct array (defined
	% only once) then length(param.totalidfields) = 1 since
	% param.totalidfields = {{id1, id2}}
	if ~strcmpi(params.ignoreTotalIDFields, 'yes')
		if ~isempty(params.totalidfields)
			if length(params.totalidfields) > 1
				totalidfields = params.totalidfields{cStruct};
			else
				totalidfields = params.totalidfields{1};
			end
			
			if ~iscell(totalidfields)
				totalidfields = { totalidfields };
			end
			
			filePos = strfind(totalidfields, 'filename');
            
            if ~iscell(filePos)
                filePos = {filePos};
            end
            
            %by Abel: if the totalidfields do not containt a fieldname
            %filename% gamble that the first totalidfields is the filename.
            if isempty(filePos) && length(params.totalidfields) > 0
                filePos = {1};
            end
			

			
			for cId = 1:length(filePos)
				if ~isempty(filePos{cId})
					fileName = retrieveField(structIn(cRow), totalidfields{cId});
					data(cRow).fileName = fileName{:};
				end
			end
		else
			data(cRow).fileName = '';
		end
	else
		data(cRow).fileName = '';
	end
	
	data(cRow).cellId = ''; % initialize; fill in in next loop
	
	%by Abel: added "cStruct" index of structure array needed to retrieve params
	data(cRow).buttonDownFcn = localGetButtonDownFcn(structIn(cRow), ...
		data(cRow).fileName, params, cStruct);
	%Print buttonDownFcn for debugging needed
	if (params.showexeceval)
		disp(data(cRow).buttonDownFcn);
	end
	
	%by Abel: prepare Z-Axes data for color plotting if needed
	if ~isempty(params.ZData)
		%Get Z-data from given leave name or directly from the users input
		if any(strfind(params.ZData{cStruct}, '$'))
			data(cRow).Z = eval(parseExpression(params.ZData{cStruct}, 'structIn(cRow)'));
		else
			data(cRow).Z = params.ZData{cStruct};
		end
		
		%now apply any previous selections (indexexpression etc,..)
		nRowsYidx = size( data(cRow).idx.Y, 1 );
		%take only idx common in all selections
		if nRowsYidx > 1
			idxY = intersect( data(cRow).idx.Y{:} );
		elseif nRowsYidx == 0
			idxY = [];
		else
			idxY = data(cRow).idx.Y{:};
		end
		% idxY is empty if no selections were made
		if ~isempty(idxY)
			data(cRow).Z = data(cRow).Z(idxY);
		end
		% set Z to empty if there is no Y-data remaining
		if isempty(data(cRow).Y)
			data(cRow).Z = [];
		end
	end
end

% derive unique cell identifier from parameter cellidfields
for cId = 1:length(params.cellidfield)
	newField = retrieveField( structIn, params.cellidfield{cId} );
	for cRow = 1:length( structIn )
		if ~ischar(newField{cRow})
			newField{cRow} = num2str(newField{cRow});
		end
		data(cRow).cellId = [data(cRow).cellId '-' newField{cRow}];
	end
end
data = data';
%%

%% localMergeRows
function data = localMergeRows(data)
% merge rows coming from the same cell; use data.cellidfields to determine
% rows coming from the same cell.
cRow = 1;
while cRow < length(data)
	% if this and next row are from the same cell, merge them
	if isequal( data(cRow).cellId, data(cRow+1).cellId )
		%merge
		data(cRow).X = [data(cRow).X data(cRow+1).X];
		data(cRow).Y = [data(cRow).Y data(cRow+1).Y];
		%delete row
		data = [data(1:cRow); data(cRow+2:end)];
	else
		% else, go on to the next row
		cRow = cRow+1;
	end
end

% then sort all vectors
for cRow=1:length(data)
	[data(cRow).X idx] = sort(data(cRow).X);
	data(cRow).Y = data(cRow).Y(idx);
end
%%

%% localGenerateAxisLabels
function label = localGenerateAxisLabels(Fields)
% Generate axes labels from the fieldnames
if  length(Fields) > 1  &&  isequal(Fields{:})
	Fields = unique(Fields);
end
NFields = length(Fields);
labelParts = VectorZip( Fields, repmat({' and '}, 1, NFields) );
labelParts(end) = [];
label = cat(2, labelParts{:});
%%

%% localGetButtonDownFcn
function buttonDownFcn = localGetButtonDownFcn(structIn, structName, params, cStruct)
%  cycle through the rows of the structure array
if ~isempty(params.execevalfnc)
	if iscell(params.execevalfnc)
		buttonDownFcn = params.execevalfnc{cStruct};
	else
		buttonDownFcn = params.execevalfnc;
	end
	
	% parse parameters enclosed by dollar signs
	dollarPos = strfind(buttonDownFcn, '$');
	if ~isequal( 0, mod( length(dollarPos), 2 ) )
		error(['Could not parse params.execevalfnc: expression shouldn''t'...
			'contain odd amount of dollar signs.']);
	end
	while dollarPos
		% first evaluate the parsed field, and convert to a string
		parsedField = retrieveField(structIn, ...
			buttonDownFcn( (dollarPos(1)+1):(dollarPos(2)-1) ));
		parsedField = parsedField{:};
		if ischar(parsedField)
			parsedField = [ '''' parsedField '''' ];
		else
			parsedField = num2str(parsedField);
		end
		% then, put in the buttonDownFcn for that row
		buttonDownFcn = [buttonDownFcn( 1:(dollarPos(1)-1) ) parsedField ...
			buttonDownFcn( (dollarPos(2)+1):end )];
		% and continue the loop...
		dollarPos = strfind(buttonDownFcn, '$');
	end
else
	if ~strcmpi(params.ignoreTotalIDFields, 'yes') && ~isempty(params.totalidfields)
		% get the information to be shown
		
		%by Abel: totalidfields are saved for each input struct separately:
		%cStruct = cell index
        if iscell(params.totalidfields{cStruct})
            totalIdFields = params.totalidfields{cStruct};
        else
            totalIdFields = params.totalidfields;
        end
        
		args = {};
		for cField = 1:length(totalIdFields)
			try
				idField = retrieveField(structIn, totalIdFields{cField});
				args = [args, {totalIdFields{cField}, idField{:} }] ;
			catch
				% do nothing
			end
		end
		% if any of the information is numeric, convert to string
		Nidx = find(cellfun('isclass', args, 'double'));
		for idx = Nidx(:)'
			args{idx} = num2str(args{idx}, '%.2f ');
		end
		
		Txt = [ sprintf(['\\bf\\fontsize{9}Datapoint from %s has following ID'...
			'parameters: \\rm'], structName), sprintf('%s : %s\n', args{:}) ];
		Txt = regexprep(Txt, '\n', ''', ''');
		buttonDownFcn = ['msgbox({''' Txt '''}'', upper(''' mfilename '''),'...
			'struct(''WindowStyle'', ''non-modal'', ''Interpreter'', ''tex''));'];
	else
		buttonDownFcn = '';
	end
end
%%

%% localSortCellData
function sortedCellData = localSortCellData(data, params)
% if cellidfields is not empty, we need to merge data from unique cells
% (draw lines between points from the same cell)
if ~isempty(params.cellidfield)
	% sort the structure on cellId
	sortedCellData = structsort(data.rows, 'cellId');
	% merge the rows coming from the same cell
	sortedCellData = localMergeRows(sortedCellData);
else
	sortedCellData = data.rows;
end
%%

%% localPlotPanel
function panel = localPlotPanel(XYPlotLines, XYPlotPoints, XFieldNames, ...
	YFieldNames, params)
% Add all generated plot objects to a KPlot panel

if isequal(params.johnson, 'yes')
	logY = 'yes';
	reverseY = 'yes';
	ylim = [1e-3 1];
	yTicks = 0 : .1 : 1;
	yTickLabels = 1 : -.1 : 0;
else
	logY = 'no';
	reverseY = 'no';
	ylim = [0 0];
	yTicks = {};
	yTickLabels = {};
end

%Apply coloring-by-ZAxes if params.ZColor is defined
% - ZAxes coloring in inplemented in Panel(), here we use a fake Panel() to
%   set up the XYPlotObjects and copy them to the final Panel() used to
%   generate the final plot.
% - Any user defined Z-Axes data is already set at this point, right after
%   the creation of the XYPlotPoints() object (see most outer loop)
ZColor = params.ZColor;
doZPlot = ~strcmpi(ZColor, 'none');
nOfPlots = length(XYPlotPoints);

if doZPlot
	myFakePanel = Panel('ZColorMap', params.ZColorMap, 'Tag', 'fakePanel', 'nodraw');
	
	for cPlot = 1:nOfPlots
		%Select what to color: params.ZColorType
		% fill: color all
		% edge: color only the edge (usefull for MarkerFaceColor = 'none')
		% face: color only face, not the edge
		if strcmpi('fill', params.ZColorType)
			XYPlotPoints(cPlot) = set(XYPlotPoints(cPlot), 'MarkerFaceColor', 'fill');
		end
		
		if strcmpi('face', params.ZColorType)
			%explicitly set MarkerEdgeColor so it won't be automatically
			%set within XYPlotObject(). The latter invoces the Z-Coloring
			%of the edge.
			XYPlotPoints(cPlot) = set(XYPlotPoints(cPlot), 'MarkerEdgeColor',...
				get(XYPlotPoints(cPlot), 'Color'));
			%MarkerFaceColor to fill invokes auto coloring within
			%XYPlotObject() and also Z-Coloring
			XYPlotPoints(cPlot) = set(XYPlotPoints(cPlot), 'MarkerFaceColor', 'fill');
		end
		
		%Set ZNaN
		XYPlotPoints(cPlot) = set(XYPlotPoints(cPlot), 'ZNaNColor', params.ZNaNColor);
		
		%Add the objects to the temporary Panel()
		myFakePanel = addPlot(myFakePanel, XYPlotPoints(cPlot), 'noredraw');
	end
	%Set Z-Scale if specified
	if ~isempty(params.ZScale)
		myFakePanel = setZScale(myFakePanel, params.ZScale);
	end
	
	%Use Panel() to create color maps. addZColors() acts on Panel() level
	%so we can calculate a global color scale for each XYPlotObject()
	myFakePanel = addZColors(myFakePanel, ZColor, 'noredraw');
	
	%Copy the XYPlotObject()'s
	XYPlotObjects = getPlotObjects(myFakePanel);
	for cPlot = 1:nOfPlots
		XYPlotPoints(cPlot) = XYPlotObjects{cPlot};
	end

end

%Now create the final Panel()
panel = Panel(...
	'logY', logY,...
	'reverseY', reverseY,...
	'ylim', ylim, ...
	'yTicks', yTicks,...
	'yTickLabels', yTickLabels,...
	'ZColor', ZColor,...
	'ZColorMap', params.ZColorMap,...
	'Tag', 'mainPanel',...
	'nodraw');
for cPlot = 1:length(XYPlotLines)
	panel = addPlot(panel, XYPlotLines(cPlot), 'noredraw');
	panel = addPlot(panel, XYPlotPoints(cPlot), 'noredraw');
end

% set labels
XLabel = localGenerateAxisLabels(XFieldNames);
YLabel = localGenerateAxisLabels(YFieldNames);
panel = set(panel, 'xlabel', XLabel, 'ylabel', YLabel, 'noredraw');
panel = set(panel, 'xlim', [params.xlim(1) params.xlim(2)], 'ylim', [params.ylim(1) params.ylim(2)], ...
	'noredraw');

% by abel: Add additional options to the panel if they were specified by
% the user
options = params.plotoptions;
for n = 1:length(options)
	if ~isempty(params.(options{n}))
		try
			panel = set(panel, options{n}, params.(options{n}), 'noredraw');
		catch errorSetPanel
			warning('SGSR:Critical', 'Error while setting plot option:%s\nLook at doc panel for valid options.\n', options{n});
			disp(getReport(errorSetPanel));
		end
	end
end

%% localPlotStats
function statsText = localPlotStats(structNames, stats)

structLabels = [];
nAnimals = [];
nCells = [];
nDots = [];
for cStruct = 1:length(structNames)
	structLabels = [structLabels structNames{cStruct} ', '];
	nCells = [nCells num2str(stats(cStruct).nCells) ', '];
	nAnimals = [nAnimals num2str(stats(cStruct).nAnimals) ', '];
	nDots = [nDots num2str(stats(cStruct).nDots) ', '];
end
structLabels = structLabels(1:end-2);
nAnimals = nAnimals(1:end-2);
nCells = nCells(1:end-2);
nDots = nDots(1:end-2);

statsText = {'\bfStatistics\rm'; ...
	['\it', structLabels, '\rm']; ...
	['\it#Animals:\rm ', nAnimals]; ...
	['\it#Cells:\rm ', nCells]; ...
	['\it#Dots:\rm ', nDots]; ...
	''; ...
	['\bfDate:\rm ', date]};


%% returnPlottedData
% by Abel: Return the plotted data in the original format of the input struct
function dataStruct = returnPlottedData (X, Y, data, xlim, ylim)
dataStruct = [];

% - Determine cleaned fieldnames
xField = cleanFields_(data.XFieldNames);
yField = cleanFields_(data.YFieldNames);
if (strcmp(xField, yField))
	yField = [yField '_modified'];
end

% - Build and return struct
for n=1:length(X)
	
	if isempty(X{n})
		[x, y] = deal(nan);
	else
		[x, y] = deal(X{n}, Y{n});
		outOfXlim = any(x < xlim(1)) || any(x > xlim(2));
		outOfYlim = any(y < ylim(1)) || any(y > ylim(2));
		if outOfXlim || outOfYlim
			[x, y] = deal(nan);
		end
	end
	
	cmd{1} = ['dataStruct(' num2str(n) ').' xField '=' '[x];'];
	cmd{2} = ['dataStruct(' num2str(n) ').' yField '=' '[y];'];
	cmd{3} = ['dataStruct(' num2str(n) ').name=[''' data.name '''];'];
	eval(cmd{1});
	eval(cmd{2});
	eval(cmd{3});
end
%%

%% cleanFields
% by Abel: remove any matlab code and delimiters from input and return clean leaf
% name.
function fieldName = cleanFields_(fieldName)
dollarPos = strfind(fieldName, '$');

if isempty(dollarPos)
	return;
end

if ~isequal( 0, mod( length(dollarPos), 2 ) )
	error('Could not parse fieldName: %s expression shouldn''t contain odd amount of dollar signs.', fieldName);
end

%remove all before first '$' and all behind second '$'
idx = [ (dollarPos(1) +1):(dollarPos(2) -1)];
fieldName  = fieldName(idx);
%%

%% fitSpline
% by Abel: fit cubic spline
%Spline-fitting requires at least two datapoints. Spline-fitting is done on all datapoints,
%but sampling is only done on intervals between included datapoints. This is to avoid
%extreme curves near the end ...
function [X, Y] =  fitSpline_(X, Y, fitRatio)

%return if X a single point
if size(X, 2) == 1
	return
end

%generate X-values for fit based on average increment of the X values and
%a ratio X data versus X fitted values
for n=1:length(X)
	avgXincrement = mean(diff(X{n}));
	minX = min(X{n});
	maxX = max(X{n});
	% 	minY = min(Y{n});
	% 	maxY = max(Y{n});
	xFit = minX:avgXincrement/fitRatio:maxX;
	
	%fit spline
	%	pp = spline(X{n}, [minY, Y{n}, maxY]);
	pp = spline(X{n}, Y{n});
	X{n} = xFit;
	Y{n} = ppval(pp,xFit);
end

%% Create gutter panels
function [ gutterPanelTop, gutterPanelBottom, gutterPanelLeft, gutterPanelRight ] = createGutterPanels( params )
gutterPanelTop = [];
gutterPanelBottom = [];
gutterPanelLeft = [];
gutterPanelRight = [];

% - Conbine gutters if more than one input structure was given
if size(params.gutterparam, 2) > 1
	gutter.top = concatCellArray_({params.gutterparam(:).top});
	gutter.bottom = concatCellArray_({params.gutterparam(:).bottom});
	gutter.left = concatCellArray_({params.gutterparam(:).left});
	gutter.right = concatCellArray_({params.gutterparam(:).right});
	%  else
	params.gutterparam = gutter;
end

mkTop = any(~cellfun(@isempty, params.gutterparam.top));
mkBottom = any(~cellfun(@isempty, params.gutterparam.bottom));
mkLeft = any(~cellfun(@isempty, params.gutterparam.left));
mkRight = any(~cellfun(@isempty, params.gutterparam.right));

xyPlotParams.Marker = { 'o'; '+' };
xyPlotParams.Color = 'k';


%!!Limits and Ticks are set within KStructPlotPage() (copied from main
%axes)
%set limits for gutter axes
gxlim = params.xlim;
gylim = params.ylim;
if gxlim(1) == -inf
	gxlim(1) = params.xmin;
end
if gxlim(2) == inf
	gxlim(2) = params.xmax;
end
if gylim(1) == -inf
	gylim(1) = params.ymin;
end
if gylim(2) == inf
	gylim(2) = params.ymax;
end

if mkTop
	xGut{1,:} = cell2mat(params.gutterparam.top{1});
	yGut{1,:} = cell2mat(params.gutterparam.top{2});
	
	%Draw + at axes ends
	xGut{2,:} = params.xlim;
	yGut{2,:} = [0 0];
		
	gutterPanelTop = Panel('Axes', false(), 'xlim', gxlim, 'Tag', 'gutterPanelTop', 'nodraw');
	gutterPanelTop = addPlot(gutterPanelTop,...
		XYPlotObject(xGut, yGut, xyPlotParams), 'noredraw');
end
if mkBottom
	xGut{1,:} = cell2mat(params.gutterparam.bottom{1});
	yGut{1,:} = cell2mat(params.gutterparam.bottom{2});
	
	%Draw + at axes ends
	xGut{2,:} = gxlim;
	yGut{2,:} = [0 0];
	
	gutterPanelBottom = Panel('Axes', false(), 'xlim', gxlim, 'Tag', 'gutterPanelBottom', 'nodraw');
	gutterPanelBottom = addPlot(gutterPanelBottom,...
		XYPlotObject(xGut, yGut, xyPlotParams), 'noredraw');
end
if mkLeft
	xGut{1,:} = cell2mat(params.gutterparam.left{1});
	yGut{1,:} = cell2mat(params.gutterparam.left{2});
	
	%Draw + at axes ends
	xGut{2,:} = [0 0];
	yGut{2,:} = gylim;
	
	gutterPanelLeft = Panel('Axes', false(), 'ylim', gylim, 'Tag', 'gutterPanelLeft', 'nodraw');
	gutterPanelLeft = addPlot(gutterPanelLeft,...
		XYPlotObject(xGut, yGut, xyPlotParams), 'noredraw');
end
if mkRight
	xGut{1,:} = cell2mat(params.gutterparam.right{1});
	yGut{1,:} = cell2mat(params.gutterparam.right{2});
	
	%Draw + at axes ends
	xGut{2,:} = [0 0];
	yGut{2,:} = gylim;
	
	gutterPanelRight = Panel('Axes', false(), 'ylim', gylim, 'Tag', 'gutterPanelRight', 'nodraw');
	gutterPanelRight = addPlot(gutterPanelRight,...
		XYPlotObject(xGut, yGut, xyPlotParams), 'noredraw');
end

%% Concat Cell Array
function combinedCell = concatCellArray_(cellArray)
combinedCell = {};
for n=1:2
	myCells = cellfun(@(x) cell2mat(x{n}), cellArray, 'uni', false);
	combinedCell{n} = {horzcat(myCells{:})};
end

%% Create edge histogram panels
function [ eHistX, eHistY ] = createHistPanels( param )
X = param.alldata.x;
Y = param.alldata.y;

gxlim = param.xlim;
gylim = param.ylim;
if gxlim(1) == -inf
	gxlim(1) = param.xmin;
end
if gxlim(2) == inf
	gxlim(2) = param.xmax;
end
if gylim(1) == -inf
	gylim(1) = param.ymin;
end
if gylim(2) == inf
	gylim(2) = param.ymax;
end

%calculate histogram
[xBins, xN, xFrac] = makeHist(X, gxlim, param.histbinnr);
[yBins, yN, yFrac] = makeHist(Y, gylim, param.histbinnr);

%prepare panels
eHistX = Panel('axes', false, 'xlim', param.xlim, 'ticksDir', 'out', 'Tag', 'eHistPanelTop', 'nodraw');
eHistY = Panel('axes', false, 'ylim', param.ylim, 'ticksDir', 'out', 'Tag', 'eHistPanelSide', 'nodraw');


%add histObjects to panels 
% # display freq in nr of elements
% % display freq in % of total elements
if strcmpi(param.histyunit, '#')
	histObj = histPlotObject(xBins, xN, 'xaxes', false);
	eHistX = addPlot(eHistX, histObj, 'noredraw');
	histObj = histPlotObject(yBins, yN, 'Dir', 'h', 'yaxes', false);
	eHistY = addPlot(eHistY, histObj, 'noredraw');
elseif strcmpi(param.histyunit, '%')
	histObj = histPlotObject(xBins, xFrac, 'xaxes', false);
	eHistX = addPlot(eHistX, histObj, 'noredraw');
	histObj = histPlotObject(yBins, yFrac, 'Dir', 'h', 'yaxes', false);
	eHistY = addPlot(eHistY, histObj, 'noredraw');	
end

%% Calc histograms
function [bins, nFreqs, nFrac] = makeHist(Data, Range, NBin)

nFrac = [];
Edges      = linspace(Range(1), Range(2), NBin+1);
BinWidth   = Edges(2)-Edges(1);
BinCenters = (Range(1)+BinWidth/2):BinWidth:(Range(2)-BinWidth/2);

[nFreqs, bins] = hist(Data, BinCenters);
nTotal = sum(nFreqs);
if (nTotal ~= 0)
    nFrac = 100*nFreqs/nTotal;
end
%%
