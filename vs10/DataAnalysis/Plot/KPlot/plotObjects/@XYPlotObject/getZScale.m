function ZScale = getZScale(XYP)
% getZScale      Returns ZScale from XYPlotObject 
%
% !ATTENTION, each XYP Object can only have 1 ZScale. If you need more,
% you'll have to add more XYP Objects to the panel.

%% ---------------- CHANGELOG -----------------------
%  Fri Aug 12 2011  Abel   
%   - Initial creation

%% Params 

%% Check Param 
if nargin > 1
    error('getZScale takes max one arguments.');
end


%% Return ZScale
ZScale = XYP.ZScale;

