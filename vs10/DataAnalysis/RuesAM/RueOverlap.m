function S = RueOverlap;
% RueOverlap - overlap analysis of all audiogram pairs
[S.CellNames,S.Nrep]=RueList;
Ncell = numel(S.CellNames);
S.isPair = nan(Ncell, Ncell); % true for physical pairs
S.Auto1 = nan(Ncell, Ncell); % autocorr of cell 1 of "pair"
S.Auto2 = nan(Ncell, Ncell); % autocorr of cell 2 of "pair"
S.AutoProduct = nan(Ncell, Ncell); % dot product of Auto1 and Auto2
S.Cross = nan(Ncell, Ncell); % Sum of Squares of crosscorr 
SFN  = fullfile('D:\Data\RueData\IBW\Audiograms\NewNames\Overlap', 'OverlapResults.mat');

for i1 = 1:Ncell,
    fn1 = S.CellNames{i1};
    AC1 = RueceptiveField(fn1);
    for i2 = i1+1:Ncell,
        fn2 = S.CellNames{i2};
        disp([fn1 '  ' fn2]);
        AC2 = RueceptiveField(fn2);
        CC = RueCross(fn1, fn2);
        isPair = isequal(lower(fn1), strrep(lower(fn2),'a','b')) ...
            || isequal(lower(fn2), strrep(lower(fn1),'a','b'));
        Auto1 = norm(AC1(:)); % RMS of autocorr of cell 1 of "pair"
        Auto2 = norm(AC2(:)); % RMS of autocorr of cell 1 of "pair"
        AutoProduct = dot(AC1(:), AC2(:)); % dot product of Auto1 and Auto2
        Cross = sum(CC(:).^2); % Sum of Squares of crosscorr
        S = local_setElements(S, i1, i2, isPair, Auto1, Auto2, AutoProduct, Cross);
        save(SFN, 'S');
    end
    pcolor(S.Cross);
    figure(gcf); pause(0.05);
end

%===========
function S = local_setElements(S, i1, i2, varargin);
for iarg=1:numel(varargin),
    FldNam = inputname(3+iarg);
    S.(FldNam)(i1,i2) = varargin{iarg};
    S.(FldNam)(i2,i1) = varargin{iarg};
end




