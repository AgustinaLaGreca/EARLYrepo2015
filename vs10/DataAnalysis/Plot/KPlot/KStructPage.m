function fig = KStructPage()
% Default page for KStructplot()
% 
% See also KStructPlot KStructPlotPage defaultPage KPlot


%% ---------------- CHANGELOG ------------------------
%  Tue Nov 22 2011  Abel   
%   - Initial creation
%  Tue Apr 24 2012  Abel   
%   - Added support for edge histograms

%% ---------------- Default parameters ---------------
defaults.name               = 'KStructplot';
%The tag 'defaultpage' is needed for the Panel() redraw function and others within Kplot
% defaults.Tag                = 'defaultpage';
%by Abel: use default resize function
defaults.ResizeFcn          = @resizefcn;

%% ---------------- Main function --------------------
fig = defaultPage(defaults.name);
set(fig, 'ResizeFcn', defaults.ResizeFcn);
end


%% ---------------- Local functions ------------------
function resizefcn(src, evt) 

%get figure children 
childHdls = get(src, 'Children');

%get main plot
mainAx = findobj(childHdls, 'Tag', 'mainPanel');
if isempty(mainAx)
	return;
end

%get main pos
mainPos = get(mainAx, 'Position');

%get gutters
statsAx = findobj(childHdls, 'Tag', 'statsPanel');
lgutAx = findobj(childHdls, 'Tag', 'gutterPanelLeft');
rgutAx = findobj(childHdls, 'Tag', 'gutterPanelRight');
bgutAx = findobj(childHdls, 'Tag', 'gutterPanelBottom');
tgutAx = findobj(childHdls, 'Tag', 'gutterPanelTop');
%resize gutters
for hdl = {lgutAx, rgutAx, bgutAx, tgutAx};
	resizeGutter_(mainPos, hdl);
end

%resize edge histograms
xHist = findobj(childHdls, 'Tag', 'eHistPanelTop');
yHist = findobj(childHdls, 'Tag', 'eHistPanelSide');
for hdl = {xHist, yHist};
	resizeEdgeHist_(mainPos, hdl);
end

%now redraw panels with RedrawOnResize=true
resizePanels_(src);
end

function resizeGutter_(mainAxPos, hdlAx)

if isempty(hdlAx)
	return;
end

gutTag = get(hdlAx{:}, 'Tag');
gutPos = get(hdlAx{:}, 'Position');

if strcmpi(gutTag, 'gutterpaneltop') || strcmpi(gutTag, 'gutterpanelbottom')
	gutPos([1,3]) = mainAxPos([1,3]);
end

if strcmpi(gutTag, 'gutterpanelleft') || strcmpi(gutTag, 'gutterpanelright')
	gutPos([2,4]) = mainAxPos([2,4]);
end

gutPos( gutPos < 0 ) = 0.01;

set(hdlAx{:}, 'Position', gutPos);
end

function resizeEdgeHist_(mainAxPos, hdlAx)

if isempty(hdlAx)
	return;
end

gutTag = get(hdlAx{:}, 'Tag');
gutPos = get(hdlAx{:}, 'Position');

if strcmpi(gutTag, 'eHistPanelTop')
	gutPos([1,3]) = mainAxPos([1,3]);
end

if strcmpi(gutTag, 'eHistPanelSide')
	gutPos([2,4]) = mainAxPos([2,4]);
end

gutPos( gutPos < 0 ) = 0.01;

set(hdlAx{:}, 'Position', gutPos);
end


function resizePanels_(src)
panels = get(src,'UserData');
for n = 1:length(panels)
	
	%by Abel: function to redraw objects within the panel() directly. This
	%was needed only for legendObjects to zoom text. Only panel() objects
	%which have their RedrawOnResize=true are redrawn.
	panels{n} = resize(panels{n});
end
set(src,'UserData',panels);
end
