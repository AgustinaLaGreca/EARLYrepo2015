function RueCorrAllPairs(dtype, Tmax);
%  function RueCorrAllPairs - compile list of Rue datafiles of given type
%   RueCorrAllPairs(dtype, Tmax); 
%     dtype is 'audiograms' (default) or 'FM'
%     Tmax is end of analysis window in ms; default 450 ms = entire rep.
%
%   Results are stored in AllCorrs subdir of datadirectory containing the
%   IBW files. 
%
%   Examples:
%       RueCorrAllPairs('audiograms', 450) % ipsi & contra treated as diff cells
%       RueCorrAllPairs('audiogramsq', 450) % ipsi & contra pooled

if nargin<1, dtype = 'audiograms'; end
if nargin<2, Tmax = 450; end

[dtype, Mess] = keywordMatch(dtype, {'audiograms' 'FM' 'audiogramsq' 'FMq'}, 'dtype');
error(Mess);

[FN, Nrep, Ncell]=ListRueData(dtype);
dtype = strrep(dtype, 'q', '');
Stats = [];
AllCorrFile = fullfile('D:\Data\RueData\IBW', dtype, 'AllCorrs','AllStats.mat');
if ~exist(AllCorrFile, 'file'), % initialize Stats
    save(AllCorrFile, 'Stats');
end
for icell_1=1:Ncell,
    fn1 = FN{icell_1};
    nrep1 = Nrep(icell_1);
    for icell_2=1:icell_1,
        fn2 = FN{icell_2};
        nrep2 = Nrep(icell_2);
        OutFileName = [fn1 '_' fn2 '_' num2str(Tmax) 'ms'];
        FFN = fullfile('D:\Data\RueData\IBW', dtype, 'AllCorrs', [OutFileName, '.txt']);
        if exist(FFN, 'file'), continue; end
        OutFileName
        CC = local_corr(fn1, fn2, nrep1, nrep2, Tmax);
        textwrite(FFN, num2str(CC), '-overwrite');
        % append Stats element. Load from file and concatenate so to be
        % robust against interruptions.
        Stats = load(AllCorrFile, 'Stats'); Stats = Stats.Stats;
        Stats = [Stats, local_stats(CC, fn1, fn2)];
        save(AllCorrFile, 'Stats');
    end
end


% ========================
function [CC] = local_corr(fn1, fn2, nrep1, nrep2, Tmax);
isAuto = isequal(fn1, fn2); % only off-diagonal correlation, irep1<irep2
CC = [];
for irep1=1:nrep1,
    sam1 = ReadRueData(fn1, irep1-1, [], Tmax);
    if isAuto, n2max = irep1-1; % off -diag only
    else, n2max = nrep2;
    end
    for irep2=1:n2max,
        sam2 = ReadRueData(fn2, irep2-1, [], Tmax);
        CC = [CC; corr(sam1, sam2)];
    end
end

function S = local_stats(CC, cell1, cell2);
S.cell_1 = cell1; S.cell_2 = cell2;
S.CC = CC;
S.mean = mean(CC); S.std = std(CC);
S.median = median(CC); S.iqr = iqr(CC);
S.tt5p_right = ttest(CC, 0, 0.05, 'right');
S.tt5p_both = ttest(CC, 0, 0.05, 'both');
S.tt1p_right = ttest(CC, 0, 0.01, 'right');
S.tt1p_both = ttest(CC, 0, 0.01, 'both');
S






