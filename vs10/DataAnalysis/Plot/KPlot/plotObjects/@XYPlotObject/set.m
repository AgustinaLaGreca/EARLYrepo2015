function P = set(P, varargin)
% SET sets a property of the XYPlotObject to a new value
%
% XYPlots have a set of properties you can retrieve or adjust. The 'get' and
% 'set' functions are used for this purpose.
%
% P = set(P, propName1, newValue1, propName2, newValue2, propName3, newValue3, ...)
% Sets the value of property propName.
%
%           P: The XYPlotObject instance.
%    propName: The property you want to get, this should be a member of the
%             'params' property of P.
%    newValue: The new value we give to the field.
%
% Properties for XYPlots: type 'help XYPlotObject'.
%
% Example:
%  Suppose P is an XYPlotObject instance, with property 
%   Color == {'r';'k'},
%  then:
%  >> get(P, 'Color')
%  ans = 
%      'r'
%      'k'
%  >> P = set(P, 'Color', {'b';'k'});
%  >> get(P, 'Color')
%  ans = 
%      'b'
%      'k'
%

%% ---------------- CHANGELOG -----------------------
%  Tue Nov 22 2011  Abel   
%  - Added check && repmat to make sure that each XYObject row (plot) has 
%    its own parameter set.  

%% Check input && repmat for the number of plots/rows if needed
nrOfPlots = nPlots(P);
%save into params 
P.params = processParams(varargin, P.params);

%check matlab fields for the correct number of rows
fNames = fieldnames(P.params.ML);
for n=1:length(fNames)
	nOfRows = size(P.params.ML.(fNames{n}), 1);
	if nOfRows < nrOfPlots
		%save as cell 
		if ~iscell(P.params.ML.(fNames{n}))
			P.params.ML.(fNames{n}) = {P.params.ML.(fNames{n})};
		end
		%repmat to nrOfPlots
		P.params.ML.(fNames{n}) = repmat(P.params.ML.(fNames{n}), nrOfPlots, 1);
	end
end

