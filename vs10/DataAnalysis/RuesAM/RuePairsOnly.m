function [M,RS]=RuePairsOnly(Dur);
% restrict cross-cell correlation to cells from pairs, contra only.

RS = RueSelectSignCells(Dur);

% select pair members, whose name start with p
ok1 = strmatch('p', {RS.CC.cell_1});
ok2 = strmatch('p', {RS.CC.cell_2});
ok  = intersect(ok1, ok2); % indices of comparisons between p* and p*
RS.CC = RS.CC(ok);

% select contralaterals, whose name start with c
nams={RS.CC.cell_1};
sman  = cellfun(@fliplr, nams , 'UniformOutput', false);
ok1 = strmatch('c', sman);
nams={RS.CC.cell_2};
sman  = cellfun(@fliplr, nams , 'UniformOutput', false);
ok2 = strmatch('c', sman);
ok  = intersect(ok1, ok2); % indices of comparisons between p* and p*
RS.CC = RS.CC(ok);

N = numel(RS.CC);
iok=[];
for ii=1:N,
    name1 = RS.CC(ii).cell_1;
    name2 = RS.CC(ii).cell_2;
    % select only comparisons of an A with a B or vice versa.
    if ismember('A', name1) && ismember('B', name2),
        iok = [iok ii];
    elseif ismember('B', name1) && ismember('A', name2),
        iok = [iok ii];
        [RS.CC(ii).cell_1, RS.CC(ii).cell_2] = swap(RS.CC(ii).cell_1, RS.CC(ii).cell_2);
    end
end
RS.CC = RS.CC(iok);
Npairs = numel(RS.CC);

RS.allNames_A = unique({RS.CC.cell_1});
RS.allNames_B = unique({RS.CC.cell_2});

M = nan(numel(RS.allNames_A), numel(RS.allNames_B));
for ii=1:Npairs,
    cc = RS.CC(ii);
    idxA = strmatch(cc.cell_1, RS.allNames_A)
    idxB = strmatch(cc.cell_2, RS.allNames_B)
    M(idxA, idxB) = cc.mean;
end
MM = [M, nan(size(M,1),1)]; MM = [MM; nan(1, size(MM,2))];
pcolor(MM);
colorbar;
xlim([-5 22])
ylim(xlim)
for ii=1:numel(RS.allNames_A),
    nam = RS.allNames_A{ii};
    text(ii+0.5, 0, nam, 'Rotation', -70 ,'VerticalAlignment', 'bo')
end
for ii=1:numel(RS.allNames_B),
    nam = RS.allNames_B{ii};
    text(0, ii+0.5, nam, 'HorizontalAlignment', 'right');
end
xplot(xlim,ylim,':', 'linewidth', 2);
title([num2str(Dur) ' ms']);







