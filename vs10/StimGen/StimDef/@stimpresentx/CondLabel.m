function L = CondLabel(P, K)
% stimPresent/CondLabel - string label for a stimulus condition
%    CondLabel(P, K) returns a label indicating the Kth stimulus
%    condition. The value and unit are taken from P.X and P.Y, that is, the
%    two varied parameters. Y may be empty, in which case it is ignored.
%
%    Note that K counts the conditions "as defined", which is
%    not necessarily the seem as the order "as played". 
%    If K is an array, CondLabel returns a cell array of strings.
%    K==0 means "all conditions."
%
%
%   See stimpresentx.

if nargin<2, K = 0; end % default: all conditions

if isequal(0,K),
    K = 1:P.Ncond;
end

if any(K>P.Ncond | K<0),
    error('Stimulus condition index is out of range.');
end

[X,Y] = deal(P.X, P.Y);
hasY = ~isempty(Y);

if numel(K)==1,
    if hasY,
        L = sprintf([X.FormatString ', ' Y.FormatString], X.PlotVal(K), Y.PlotVal(K));
    else,
        L = sprintf(X.FormatString, X.PlotVal(K));
    end
else, % insert tabs and use them to separate 
    if hasY,
        L = sprintf([X.FormatString ', ' Y.FormatString '\t'], VectorZip(X.PlotVal(K), Y.PlotVal(K)) );
    else,
        L = sprintf([X.FormatString '\t'], X.PlotVal(K));
    end
    L = Words2cell(L, char(9));
end








