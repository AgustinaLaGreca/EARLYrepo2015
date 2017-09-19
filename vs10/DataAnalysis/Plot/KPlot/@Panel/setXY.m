function myPanel = setXY(myPanel, varargin)
% SET sets a property of All XYObjects of the Panel instance to a new value
%
% Panels have a set of properties you can retrieve or adjust. The 'get'
% and 'set' functions are used for this purpose.
%
% panel = set(panel, propName1, newValue1, propName2, newValue2, propName3, newValue3, ...)
% Sets the value of property propName.
%
%       panel: The Panel instance.
%    propName: The property you want to set.
%    newValue: The new value we give to the property.
%
% Properties for Panel: see 'help Panel'.
%
% Example:
%  Suppose panel is a Panel instance, with property Title == 'Figure 1',
%  then:
%  >> get(panel, 'Title')
%  ans =
%      'Figure 1'
%  >> panel = set(panel, 'Title', 'Experiment Results');
%  >> get(panel, 'Title')
%  ans =
%      'Experiment Results'

%% Check Params 
% last option redraw?
noredraw = 0;
if isequal('noredraw', lower(varargin{end}))
    noredraw = 1;
	varargin = {varargin{1:end-1}};
end

%Get index of the XYObjects we want to apply the property change to
nrOfPlots = nPlots(myPanel);
if strcmpi(varargin{1}, 'all')
	plotIndex = 1:nrOfPlots;
	varargin = {varargin{1:end-1}};
elseif isnumeric(varargin{1})
	plotIndex = varargin{1};
	varargin = {varargin{1:end-1}};
else
	plotIndex = 1:nrOfPlots;
end

%% Set Values 
% params = processParams(varargin, myPanel.params);
%loop over XYObjects
for n=plotIndex
	myPanel.plotObjects{n} = set(myPanel.plotObjects{n}, varargin{:});
end


if ~noredraw
    myPanel = redraw(myPanel);
end