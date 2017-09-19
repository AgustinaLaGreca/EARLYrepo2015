function HistPlotObject = redraw(hdl, HistPlotObject, varargin)
% REDRAW Redraws the given plot on the given handle.
%
% HistPlotObject P = Redraw(HistPlotObject P, dbl hdl[, int startColor])
% The plot P is being plotted on the given figure or axes.
%
%           P: The HistPlotObject which is being plotted.
%         hdl: The handle where HistPlotObject is being plotted (typically a figure
%              or an axes).

% Created by: Kevin Spiritus
% Last edited: December 4th, 2006

%% ---------------- CHANGELOG -----------------------
%  Tue Apr 24 2012  Abel: rewrite   
%   - add support for more params/properties
%   - add support for barh()


%% parameters
if isequal(2, nargin)
    startColor = 1;
end

if nargin > 3
    error('Too many arguments.');
end

%% check given handle
try
    findobj(hdl);
catch
    error('Given handle for drawing HistPlotObject is invalid. Type ''help redraw'' for more information.');
end

%% Then plot
axes(hdl);
% bar(hdl, HistPlotObject.XData, HistPlotObject.YData, HistPlotObject.params.Color);
plotStr = 'HistPlotObject.plotHdl = ';
if strcmpi(HistPlotObject.params.Dir, 'h')
	plotStr = [plotStr, 'barh(hdl, '];
else
	plotStr = [plotStr, 'bar(hdl, '];
end
X = HistPlotObject.XData;
Y = HistPlotObject.YData;

color = HistPlotObject.params.Color;
if isnumeric(color)
    plotStr = [plotStr, 'X, Y,''FaceColor''' , ',[' , num2str(color), ']', ', ' ];
else
    plotStr = [plotStr, 'X, Y, ', '''', color, '''', ', '];
end

fNames = fieldnames(HistPlotObject.params.ML);
for n = 1:length(fNames)
	plotStr = [ plotStr, sprintf('''%s'', ', fNames{n}, HistPlotObject.params.ML.(fNames{n})) ];
end
plotStr = plotStr(1:end-2);
plotStr = [plotStr, ');'];

%save axes param
axParam = get(hdl);

%eval 
eval(plotStr);

%% Reset parent poperties which were overwritten by bar() ?? 
for n = 1:length(HistPlotObject.params.inherit)
	set(hdl, HistPlotObject.params.inherit{n}, axParam.(HistPlotObject.params.inherit{n}));
end

if ~HistPlotObject.params.xaxes
	set(hdl, 'XTick', []);
end
if ~HistPlotObject.params.yaxes
	set(hdl, 'YTick', []);
end

	


