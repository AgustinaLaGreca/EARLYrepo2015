function newHistPlot = HistPlotObject(XData, YData, varargin)
% HistPLOTOBJECT Creates a HistPlotObject instance.
%
% An HistPlotObject contains data and information for a histogram plot.
% Plotting happens using the redraw method. After construction the
% HistPlotObject can be added to a Panel object. 
%
% newHistPlot = HistPlotObject(XData, YData, paramList) 
% Creates an HistPlotObject object from the given vectors. Parameter pairs
% can be given in paramList.
%
%       XData: XData for the histogram. It is a row vector.
%       YData: YData for the histogram. It is a row vector.
%       Color: Color for the histogram; is be a single value

% Created by: Kevin Spiritus
% Last adjusted: April 30th, 2007


%% ---------------- CHANGELOG -----------------------
%  Tue Apr 24 2012  Abel   
%   - add support for more ML params/properties

%% Parameter checking
if nargin < 2
    error('HistPlotObject expects at least two arguments. Type ''help HistPlotObject'' for more information.');
end

if ~isnumeric(XData) || ~isequal(1, size(XData, 1))
    error('Wrong format of the arguments. Type ''help HistPlotObject'' for more information.');
end

if ~isnumeric(YData) || ~isequal(1, size(YData, 1))
    error('Wrong format of the arguments. Type ''help HistPlotObject'' for more information.');
end

paramsIn = varargin; % these will be checked later, by processParams

%% Assign data
newHistPlot.XData = XData;
newHistPlot.YData = YData;
newHistPlot.plotHdl = [];

%% Handle parameters
% define structure of parameters: just fill with defaults
% newHistPlot.params.Color = 'b';

%by Abel: Matlab properties
newHistPlot.params.ML.Visible = 'on';

%by Abel: internal properties 
newHistPlot.params.Dir = 'v';
newHistPlot.params.Color = 'b';
newHistPlot.params.xaxes = true;
newHistPlot.params.yaxes = true;
newHistPlot.params.inherit = {...
	'XDir', 'Visible', 'Tag'};




% Parameters given as arguments are placed in the params structure.
% Other entries remain 'def'. This def value is replaced when plotting.
newHistPlot.params = processParams(paramsIn, newHistPlot.params);

%% Return object
newHistPlot = class(newHistPlot, 'histPlotObject');