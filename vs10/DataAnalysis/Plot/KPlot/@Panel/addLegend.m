function panel = addLegend(panel, legend, doRedraw)
% addLegend adds a legend to the given Panel, and redraws.
%
% panel = addLegend(Panel panel, legendObject legend, 'noredraw') 
% Adds legend to panel.
%     
%      panel : A Panel instance we want to add a legend object to
%    legend : The legend object we want to add to panel
%  'noredraw': (optional) If this string is added to the argument list, the
%              plot is not redrawn. This might speed up things when adding
%              much information to a Panel.
%
% Returns: panel, with legend added to it.
% 
% Example: 
%  >> legend = legendObject('konijn');
%  >> panel = Panel;
%  >> panel = addLegend(panel, legend);

% Created by: Abel Jonckheer



%% ---------------- CHANGELOG -----------------------
%  Mon Nov 21 2011  Abel   
%   Bugfix, the UserData field of defaulPage() was not updated with the
%   current Panel() instances if doRedraw was set.

%% ---------------- TODO -----------------------
% Update documentation 




%% Check params
switch nargin
    case 2
        doRedraw = 1;
    case 3
        if isequal('noredraw', doRedraw)
            doRedraw = 0;
        else
            error('Only two argumens expected, unless an extra ''noredraw'' is given. Type ''help addLegend'' for more information.');
        end
    otherwise
        error('This function expects two or three arguments. Type ''help addLegend'' for more information.');
end

%% add the plot and redraw
if ~( isequal('legendObject', class(legend)) || isequal('HeaderObject', class(legend)) )
    error('addLegend only accepts legendObjects (or HeaderObjects).');
end

panel.legendObjects{end+1} = legend;
if isequal(1,doRedraw)
	%Special case: If the panel was already drawn (and registered), the redraw done here
	%will not add the panel to the UserData field of the page. This may
	%lead to problems with the resize() function in the defaultPage().
	%We'll have to re-add the panel again.
	wasReg = panel.isRegistered;
	%draw panel
    panel = redraw(panel);
	%update userdata field
	if wasReg
		[userdata, page] = getUserData_();
		if ~isempty(userdata)
			%get panels in userdata
			for n=1:length(userdata)
				panelHdl = getHdl(userdata{n});
				%replace panel in userdata with same hdl
				if panelHdl == panel.hdl;
					userdata{n} = panel;
				end
			end
			%reset panels in userdata
			set(page,'UserData',userdata);
		end
	end
end

function page = getPage_()
page = [];
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
end
function [userdata, page] = getUserData_()
userdata = {};
page = getPage_();
if ~isempty(page)
	userdata = get(page,'UserData');
end