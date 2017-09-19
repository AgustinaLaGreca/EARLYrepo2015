function CatPlot(T, M, I, flag);
% CatPlot - plot categorized columns of matrix and derived templates
%   CatPlot(T, M, I), where T is a X-axis array, M is a matrix, 
%   and I is a cell array containing index arrays that select columns of M,
%   produces two graphs in a single figure. The first graph plots all 
%   columns of M against vector T, while using the same color for all 
%   columns of the same category. The second graph plots the corresponding 
%   category templates. T must have length equal to size(M,1) or [], in 
%   which case 1:size(M,1) is used.
%
%   Note:
%   I is typically the output of Categorize.
%
%   Example 
%   M = randn(10, 31)
%   Crit = std(M); % categorization of M's columns based on their STD 
%   I = categorize(Crit, 6); % six ~equal-sized categories of increasing STD
%   CatPlot([], M, I)
%    
%
%   See also Categorize, CatPlot, MEDIAN.

if nargin<4, flag=''; end

if isempty(T),
    T = 1:size(M,1);
end
T = T(:);
if ~isequal(numel(T), size(M,1)),
    error('Length of T must be equal to # rows of M.');
end

set(gcf,'units', 'normalized', 'position', [0.198 0.42 0.773 0.463])
% templates
hax2 = subplot(1,2,2);
[TM, Res, Ndist] = MultiMedian(M,I, 1);
[dum, isort] = sort(Ndist);
%f7; hist(Ndist,33); pause; delete(7);
plot(T, TM, 'linewidth', 2);
for ii=1:numel(I),
    LegStr{ii} = ['N=' num2str(numel(I{ii}))];
end
YL = ylim;
legend(LegStr);
% indiv waves
hax1 = subplot(1,2,1);
cla;
if isequal('residue', lower(flag)),
    M = Res;
end
for icat=numel(I):-1:1,
    Icur = I{icat};
    if ~isempty(Icur),
        xplot(T, M(:, Icur), lico(icat));
    end
end
xplot(T, median(M,2), 'linewidth',3, 'color', 0.8*[1 1 1]);
% xplot(T, M(:, isort(end-3:end)), 'linewidth', 3)
linkaxes([hax1 hax2], 'y');
ylim(YL);

set(figure,'units', 'normalized', 'position', [0.701 0.0525 0.266 0.346]);
hist(Ndist, 33)



