function L = CondLabel(P, K, FN)
% stimPresent/CondLabel - string label for a stimulus condition
%    CondLabel(P, K) returns a label indicating the Kth stimulus
%    condition. The value and unit are taken from P.X.
%    Note that K counts the conditions "as defined", which is
%    not necessarily the seem as the order "as played". 
%    If K is an array, CondLabel returns a cell array of strings.
%    K==0 means "all conditions."
%
%    CondLabel(P, K, 'Y') uses independent variable Y instead of the
%    default X, etc.
%
%   See stimpresent.

if nargin<2, K = 0; end % default: all conditions
if nargin<3, FN = 'X'; end % default: indep var X

X = P.(FN);
if isequal(0,K),
    K = 1:numel(X.PlotVal);
end

if numel(K)==1,
    L = sprintf(X.FormatString, X.PlotVal(K));
else, % insert tabs and use them to separate 
    L = sprintf([X.FormatString '\t'], X.PlotVal(K));
    L = Words2cell(L, char(9));
end








