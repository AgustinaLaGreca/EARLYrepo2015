function [I, N, Bounds, LegStr] = categorize(X, Bounds, LegVarname);
% categorize - divide array elements over bins
%   I = categorize(X, B) returns cell array I whose elements contain
%   indices that divide the elements of X into categories specified by
%   the boundaries in B:
%     B(n) < X(I{k}(n)) < B(n+1)
%   Use -inf and inf to get single-sided unbounded categories.
%   By convention, NaN-valued elements of X fall outside any category.
%
%   If B is a single positive number, the X values are sorted into B 
%   categories having approximately equal number of members. If B is a
%   single negative number, a number of ~equal-size categories is chosen
%   such that each contains approximately abs(B) members.
%
%   [I, N, B] = categorize(X, B) also returns the number N of elements in 
%   each category and the true Bounds values.
%
%   [I, N, B, L] = categorize(X, B, 'Foo') also returns cell string L that
%   describes each categories, e.g.,   '1.2 < Foo < 8 (N=42)'
%
%   See also Betwixt, MultiMedian, HIST.

if nargin<3, LegVarname = 'X'; end
LegVarname = cellify(LegVarname);
if numel(LegVarname)>1,
    LegVarUnit = [' ' LegVarname{2}];
    LegVarname = LegVarname{1};
else,
    LegVarUnit = '';
    LegVarname = LegVarname{1};
end

% remove nans, but remember were they were
[Xclean, inonnan] = denan(X);
% add a little noise to Xclean to eliminate coincident X values. Noise is small
% enough not to change true (nonzero) differences between X values.
[Xs, isort] = sort(Xclean);
dx = abs(diff(Xs));
mindx = min(dx(dx>0));
Xclean = Xclean + 0.5*mindx*rand(size(Xclean));
Xclean = sort(Xclean);

Ncat = numel(Bounds)-1;
if isequal(0, Ncat), % divide into ~equal-size categories (see help text)
    if Bounds<0, 
        Bounds = round(numel(Xclean)/abs(Bounds));
    end
    [Xs, isort] = sort(Xclean);
    Lx = numel(Xclean);
    Ncat = Bounds;
    x = linspace(0,1,Lx);
    nx = linspace(0,1,Ncat+1);
    Bounds = interp1(x, Xs, nx);
    Bounds([1 end]) = [-inf inf];
end
I = cell(1,Ncat); N = zeros(1,Ncat); LegStr = {};
for icat=1:Ncat,
    idx = find(betwixt(Xclean,Bounds(icat:icat+1)));% idx is Xclean-indexinf array ...
    I{icat} = inonnan(idx);% this is the X-indexing variant
    N(icat) = numel(I{icat});
    LegStr{icat} = [num2str(Bounds(icat),2) LegVarUnit ' < ' LegVarname ' < ' num2str(Bounds(icat+1),2) LegVarUnit];
end







