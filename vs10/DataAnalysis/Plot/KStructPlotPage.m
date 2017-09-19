function [ axH, param ] = KStructPlotPage(KStructplotparam, param, doResize)
%KSTRUCTPLOTPAGE		Interface between KStructplot() and PlotPage() 
%

%% ---------------- CHANGELOG ------------------------
%  Tue Nov 22 2011  Abel   
%   - Initial creation
%  Tue Apr 24 2012  Abel   
%   - Added support for edge histograms 


%% ---------------- Default parameters ---------------
axH = 0;

%Set some defaults for PlotPage
defparam = plotPage();
defparam.statsstring = '';
param = updatestruct(defparam, param);
param.nrofrows = 12;
param.nrofcols = 12;
param.pageFunction = 'KStructPage()';
param.plotpositions = {};
param.panelObjects = [];
% param.intraHmarginFactor = 0.98;
param.dateStringObject = dateStringPanel('text', ...
	[datestr(now), '     <', KStructplotparam.inputstring, '>     '],...
	'RedrawOnResize', true, 'Tag', 'dateStringPanel' );  

%Set options for KStructPlotPage
addStats = ~isempty(param.statsstring);
addTopGutter = ~isempty(param.gutterpaneltop);
addBottomGutter = ~isempty(param.gutterpanelbottom);
addLeftGutter = ~isempty(param.gutterpanelleft);
addRightGutter = ~isempty(param.gutterpanelright);
addXeHist = ~isempty(param.ehistx);
addYeHist = ~isempty(param.ehisty);
panelSequence = {};

%Used in page resize function
if nargin < 3
	doResize = false;
end

%% ---------------- Main function --------------------
%If doResize, just redraw the gutters and return 
if doResize
	%resize gutter
	resizeGutter_(param,addRightGutter,addLeftGutter,addTopGutter,addBottomGutter);
	%resize edge histogram
	resizeEdgeHist_(param,addXeHist,addYeHist);
	return;
end

%Reset Position 
param.mainpanel = set(param.mainpanel, 'Position', [0 0 0 0], 'noredraw');

%Main positions
mainPos = zeros(param.nrofrows, param.nrofcols);
for nRow = 1:param.nrofrows;
	mainPos(nRow, :) = 1+ (nRow -1) * param.nrofcols : (nRow) * param.nrofcols;
end

%Add statistics
if addStats
	% statsText && Panel
	%put the textbox NorthWest to align it with the left side of the
	%Panel() axes
	statsBox = textBoxObject(param.statsstring, 'Position', 'NorthWest');
	statsPanel = Panel('axes', false, 'Tag', 'statsPanel', 'nodraw');
	statsPanel = addTextBox(statsPanel, statsBox, 'noredraw');	
	%get stat positions -> take last colum
	statsPos = torow(mainPos(:, param.nrofcols));
	%remove 1 col in mainpanel positions
	mainPos = mainPos(:, 1:param.nrofcols -1);
	%save for PlotPage()param.nrofcols
	param.plotpositions = { statsPos };
	param.panelObjects = [ statsPanel ];
	panelSequence{1, end+1}= 'statspanel';
end

%Add Gutters
if addRightGutter
	%take most right col 
	rGutPos = torow(mainPos(:,end));
	%remove 1 col in mainpanel positions
	mainPos = mainPos(:, 1:end -1);
	%set tag for KStructPage() resize function
	param.gutterpanelright = set(param.gutterpanelright, 'Tag', 'gutterPanelRight','noredraw');
	%set hidden for now -> reset after redraw
	param.gutterpanelright = setXY(param.gutterpanelright, 'Visible', 'off', 'noredraw');
	%save for PlotPage()param.nrofcols
	param.panelObjects = [ param.panelObjects param.gutterpanelright ];
	param.plotpositions{1, end+1} = torow(rGutPos);
	panelSequence{1, end+1}= 'gutterpanelright';
end
if addLeftGutter
	%take most left col 
	lGutPos = torow(mainPos(:,1));
	%remove 1 col in mainpanel positions
	mainPos = mainPos(:, 2:end);
	%set tag for KStructPage() resize function
	param.gutterpanelleft = set(param.gutterpanelleft, 'Tag', 'gutterPanelLeft','noredraw');
	%set hidden for now -> reset after redraw
	param.gutterpanelleft = setXY(param.gutterpanelleft, 'Visible', 'off', 'noredraw');
	%save for PlotPage()param.nrofcols
	param.panelObjects = [ param.panelObjects param.gutterpanelleft ];
	param.plotpositions{1, end+1} = torow(lGutPos);
	panelSequence{1, end+1}= 'gutterpanelleft';
end
if addTopGutter
	%take most top row 
	tGutPos = torow(mainPos(1,:));
	%remove 1 row in mainpanel positions
	mainPos = mainPos(2:end,:);
	%set tag for KStructPage() resize function
	param.gutterpaneltop = set(param.gutterpaneltop, 'Tag', 'gutterPanelTop','noredraw');
	%set hidden for now -> reset after redraw
	param.gutterpaneltop = setXY(param.gutterpaneltop, 'Visible', 'off', 'noredraw');
	%save for PlotPage()param.nrofcols
	param.panelObjects = [ param.panelObjects param.gutterpaneltop ];
	param.plotpositions{1, end+1} = torow(tGutPos);
	panelSequence{1, end+1}= 'gutterpaneltop';
end	
if addBottomGutter
	%take most bottom row 
	bGutPos = torow(mainPos(end,:));
	%remove 1 col in mainpanel positions
	mainPos = mainPos(1:end-1,:);
	%save for PlotPage()param.nrofcols
	%set tag for KStructPage() resize function
	param.gutterpanelbottom = set(param.gutterpanelbottom, 'Tag', 'gutterPanelBottom','noredraw');
	%set hidden for now -> reset after redraw
	param.gutterpanelbottom = setXY(param.gutterpanelbottom, 'Visible', 'off', 'noredraw');
	param.panelObjects = [ param.panelObjects param.gutterpanelbottom ];
	param.plotpositions{1, end+1} = torow(bGutPos);
	panelSequence{1, end+1} = 'gutterpanelbottom';
end	

%Add edge histograms
if addXeHist
	%take top row
	xHistPos = torow(torow(mainPos(1,:)));
	%remove top row from main
	mainPos = mainPos(2:end,:);
	%set tag for KStructPage() resize function
	param.ehistx = set(param.ehistx, 'Tag', 'eHistPanelTop', 'noredraw');
	param.panelObjects = [ param.panelObjects param.ehistx];
	param.plotpositions{1, end+1} = torow(xHistPos);
	panelSequence{1, end+1} = 'ehistpaneltop';
end
if addYeHist
	%take left col
	yHistPos = torow(mainPos(:,end));
	%rm left col from main
	mainPos = mainPos(:, 1:end -1);	
	%set tag for KStructPage() resize function
	param.ehisty = set(param.ehisty, 'Tag', 'eHistPanelSide', 'noredraw');
	param.panelObjects = [ param.panelObjects param.ehisty];
	param.plotpositions{1, end+1} = torow(yHistPos);
	panelSequence{1, end+1} = 'ehistpanelleft';
end

%Add main panel to page
param.panelObjects = [ param.panelObjects param.mainpanel ];
param.plotpositions{1, end+1} = torow(mainPos);
panelSequence{1, end+1}= 'mainpanel';

%Now plot
[ axH, param ] = plotPage(param);

%Reset the panels after plotting 
for nObj = 1:length(panelSequence)
	param.(panelSequence{nObj}) = param.panelObjects(nObj);
end

%Resize gutter
resizeGutter_(param,addRightGutter,addLeftGutter,addTopGutter,addBottomGutter);

%Resize eHist
resizeEdgeHist_(param, addXeHist, addYeHist);

%Add legend
legendHdls = getFirstHandles(param.mainpanel);
% keep even entries only (they contain the markers)
legendHdls = legendHdls(2:2:end);
[legendHdls, nanIdx] = denan(legendHdls);
allHdls = getHandles(param.mainpanel);
labels = repmat({''}, 1, length(allHdls));
idx = ismember(allHdls, legendHdls);

%Needed for empty legend param 
if ~iscell(KStructplotparam.legend)
	KStructplotparam.legend = {KStructplotparam.legend};
end

if all(cellfun(@isempty, KStructplotparam.legend))
	labels(idx) = KStructplotparam.structnames(nanIdx);
else
	labels(idx) = KStructplotparam.legend;
end

if isempty(KStructplotparam.legendParam)
	lObj = legendObject('textlabels', labels);
else
	lObj = legendObject('textlabels', labels, KStructplotparam.legendParam{:});
end

%Set MainAxes to be active 
axes(axH(2));
%add legend and redraw
param.mainpanel = addLegend(param.mainpanel, lObj);

end

%% ---------------- Local functions ------------------
function resizeGutter_(param,addRightGutter,addLeftGutter,addTopGutter,addBottomGutter)
	%Set positions of gutter
	mainAx = getHdl(param.mainpanel);
	mainPos = get(mainAx, 'Position');
 	mainXtick = get(mainAx, 'XTick');
	mainYtick = get(mainAx, 'YTick');
	mainXlim = get(mainAx, 'XLim');
	mainYlim = get(mainAx, 'YLim');
	
	if addRightGutter
        try % TODO: remove
            gHdl = getHdl(param.gutterpanelright);
            myPos = get(gHdl, 'Position');
            set(gHdl, 'Position', [myPos(1), mainPos(2), myPos(3), mainPos(4)]);
            set(gHdl, 'YLim', mainYlim);
            set(gHdl, 'YTick', mainYtick);
            %set all line handles to 'on'
            set(getHandles(param.gutterpanelright), 'Visible', 'on');
        end
	end
	if addLeftGutter
		gHdl = getHdl(param.gutterpanelleft);
		myPos = get(gHdl, 'Position');
		set(gHdl, 'Position', [myPos(1), mainPos(2), myPos(3), mainPos(4)]);
		set(gHdl, 'YLim', mainYlim);
		set(gHdl, 'YTick', mainYtick);
		%set all line handles to 'on'
		set(getHandles(param.gutterpanelleft), 'Visible', 'on');
	end
	if addTopGutter
		gHdl = getHdl(param.gutterpaneltop);
		myPos = get(gHdl, 'Position');
		set(gHdl, 'Position', [mainPos(1), myPos(2), mainPos(3), myPos(4)]);
		set(gHdl, 'XLim', mainXlim);
		set(gHdl, 'XTick', mainXtick);
		%set all line handles to 'on'
		set(getHandles(param.gutterpaneltop), 'Visible', 'on');
	end
	if addBottomGutter
		gHdl = getHdl(param.gutterpanelbottom);
		myPos = get(gHdl, 'Position');
		set(gHdl, 'Position', [mainPos(1), myPos(2), mainPos(3), myPos(4)]);
		set(gHdl, 'XLim', mainXlim);
		set(gHdl, 'XTick', mainXtick);
		%set all line handles to 'on'
		set(getHandles(param.gutterpanelbottom), 'Visible', 'on');
	end
end
function resizeEdgeHist_(param,addXeHist,addYeHist)
	%Set positions of egde histogram
	mainAx = getHdl(param.mainpanel);
	mainPos = get(mainAx, 'Position');
 	mainXtick = get(mainAx, 'XTick');
	mainYtick = get(mainAx, 'YTick');
	mainXlim = get(mainAx, 'XLim');
	mainYlim = get(mainAx, 'YLim');
	
	if addXeHist
		hHdl = getHdl(param.ehistpaneltop);
		myPos = get(hHdl, 'Position');
		set(hHdl, 'Position', [mainPos(1), myPos(2), mainPos(3), myPos(4)]);
		set(hHdl, 'XLim', mainXlim);
		set(hHdl, 'XTick', []);
		%set all handles to 'on'
		set(getHandles(param.ehistpaneltop), 'Visible', 'on');
	end
	
	if addYeHist
		hHdl = getHdl(param.ehistpanelleft);
		myPos = get(hHdl, 'Position');
		set(hHdl, 'Position', [myPos(1), mainPos(2), myPos(3), mainPos(4)]);
		set(hHdl, 'YLim', mainYlim);
		set(hHdl, 'YTick', []);
		%set all handles to 'on'
		set(getHandles(param.ehistpanelleft), 'Visible', 'on');
	end
end
