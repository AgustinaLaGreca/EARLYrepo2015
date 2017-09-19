function sscatter(X, Y, M, Bounds, Mrk);
% sscatter - slow scatter plot
%   qscatter(X, Y, M, mB, mS)
%    X: quantity at X-axis
%    Y: field name of quantity at Y-axis
%    M: quantitity determining marker colors
%    mB: array holding marker boundaries. If mB is omitted, but M is
%        specified, omitted, each value of S.Fm results in a different
%        marker.
%    mS: struct array holding marker props corresponding to subgroups
%        determined by Fm and mB.
%
%    Omitting an input arg has the same effect as specifying [].
%    qscatter works both for struct arrays and single structs whos fields
%    are arrays.

if nargin<3, M=[]; end % default: same marker colors
if nargin<4, Bounds=[]; end % default: same marker colors
if nargin<5, Mrk=[]; end % default: same marker colors

N = numel(X);

LegStr = {};
if isempty(M),
    Ngroup = 1;
    SubDiv = {1:N}; % all points in single group
else,
    allMval = unique(M);
    if isempty(Bounds),
        Ngroup = numel(allMval);
        for igroup=1:Ngroup,
            m = allMval(igroup); % the value of S.Fm for this group
            SubDiv{igroup} = find(M==m);
            if isnumeric(M) || islogical(M),
                LegStr = [LegStr, [Fm '=' num2str(m)]];
            elseif ischar(M);
                LegStr = [LegStr, [Fm '=' m]];
            else,
                LegStr = [LegStr, ['Group # ' num2str(igroup)]];
            end
        end
    else,
        [SubDiv, dum, Bounds, LegStr] = categorize(M, Bounds, 'X');
        Ngroup = numel(Bounds)-1;
    end
end

if isempty(Mrk),
    for ii=1:Ngroup,
        Mrk(ii).color = getfield(lico(ii),'color');
        Mrk(ii).marker = ploma(ii);
        Mrk(ii).linestyle = 'none';
    end
end

clf;

for igroup=1:Ngroup,
    idx = SubDiv{igroup};
    h = xplot(X(idx),Y(idx), Mrk(igroup));
    %IDpoints(h, S, idx, @(x,y)['event #'  num2str(y)], 'view trace', @plotSpikeMetrics);
end

if ~isempty(LegStr),
    legend(LegStr);
end
%xlabel(Fx); ylabel(Fy);

iok =  ~isnan(X) & ~isnan(Y);
RR = corrcoef(X(iok),Y(iok)); RR=RR(1,2);
title(['R = ' num2str(RR,2)]);









