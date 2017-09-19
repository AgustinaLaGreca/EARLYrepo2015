function XYP = setZScale(XYP, zScale)
%setZScale       Sets ZScale for each plot in plotIndex (defaults to all). 
%
% !ATTENTION, each XYP Object can only have 1 ZScale. If you need more,
% you'll have to add more XYP Objects to the panel.

%% ---------------- CHANGELOG -----------------------
%  Fri Aug 12 2011  Abel   
%   - Initial creation 
%  Tue Nov 22 2011  Abel   
%   - Removed plotIndex input since each XYObject can only contain a single
%     ZScale + added format checking 

%% Check Param 
if nargin < 2
    error('setZScale takes 2 arguments.');
end

%format checking
% - must be double
if iscell(zScale)
	zScale = zScale{:};
end
% - must be single row
if size(zScale, 1) > 1
	error('ZScale must be a single row vector');
end

%% Set ZScale
XYP.ZScale = zScale;
