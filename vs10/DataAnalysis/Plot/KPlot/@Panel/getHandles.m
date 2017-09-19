function hdls = getHandles(panel, getWhat)
% GETHANDLES returns the handles of the Panel object
%
% hdls = getHandles(Panel panel)
% Returns an array with all plot handles, in the order they were added to
% the panel. If a plot object in panel contains multiple plots, they will
% be returned in the order they were passed to the constructor of that plot
% object.

% Created by: Kevin Spiritus
% Last edited: April 26th, 2007

%% ---------------- CHANGELOG -----------------------
%  Tue Jan 24 2012  Abel
%   - added 'getWhat' option

hdls = [];
if nargin == 1
	myObjects = panel.plotObjects;
elseif strcmpi(getWhat, 'text')
	myObjects = panel.textObjects;
elseif strcmpi(getWhat, 'legend')
	myObjects = panel.legendObjects;
elseif strcmpi(getWhat, 'all')
	myObjects = [panel.plotObjects, panel.textObjects, panel.legendObjects];
end

hdls = local_loopOverObjects(myObjects);

end

function hdls = local_loopOverObjects(objects)

hdls = [];
nObjects = size(objects, 2);
for i = 1:nObjects
	hdls = [hdls getHandles(objects{i})];
end
end



