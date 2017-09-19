function qscatter(S, Fx, Fy, Fm, Bounds, Mrk);
% qscatter - quick scatter plot
%   qscatter(S, Fx, Fy, Fm, mB)
%     S: struct
%    Fx: field name of quantity at X-axis
%    Fy: field name of quantity at Y-axis
%    Fm: field name of quantitity determining marker colors
%    mB: array holding marker boundaries. If mB is omitted, but Fm is
%        specified, omitted, each value of S.Fm results in a different
%        marker.
%    mS: struct array holding marker props corresponding to subgroups
%        determined by Fm and mB.
%
%    Omitting an input arg has the same effect as specifying [].
%    qscatter works both for struct arrays and single structs whos fields
%    are arrays.

if nargin<4, Fm=[]; end % default: same marker colors
if nargin<5, Bounds=[]; end % default: same marker colors
if nargin<6, Mrk=[]; end % default: same marker colors

% get data
X = [S.(Fx)];
Y = [S.(Fy)];
N = numel(X);

LegStr = {};
if isempty(Fm),
    Ngroup = 1;
    SubDiv = {1:N}; % all points in single group
else,
    M = [S.(Fm)];
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
        Bounds = [-inf Bounds(:)' inf];
        Ngroup = numel(Bounds)-1;
        for igroup=1:Ngroup,
            lobound = Bounds(igroup);
            hibound = Bounds(igroup+1);
            SubDiv{igroup} = find(betwixt(M, lobound, hibound));
            if ~isempty(SubDiv{igroup}),
                LegStr = [LegStr, [num2str(lobound) '<' Fm '<' num2str(hibound)]];
            end
        end
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
    IDpoints(h, S, idx, @(x,y)['event #'  num2str(y)], 'view trace', @plotSpikeMetrics);
end
if ~isempty(LegStr),
    legend(LegStr);
end
xlabel(Fx);
ylabel(Fy);

iok =  ~isnan(X) & ~isnan(Y);
RR = corrcoef(X(iok),Y(iok)); RR=RR(1,2);
if isfield(S, 'ExpID'),
    IDstr = [S.ExpID '/' S.RecID '/' num2str(S.icond) ' ---- ' num2str(S.Xval) ' ' S.Xunit ' --- '];
else,
    IDstr = '';
end
title([IDstr 'R = ' num2str(RR,2)]);









