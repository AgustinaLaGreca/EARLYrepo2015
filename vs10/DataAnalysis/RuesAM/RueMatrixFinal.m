function Stats = RueMatrixFinal(kw)
% RueMatrixFinal - final crosscor analysis of Rue's audiogram data

switch(kw),
    case 'pooled',
        % load I&C-pooled crosscor data
        load('D:\Data\RueData\IBW\Audiograms\Corr_450ms_pooled\AllStats_450ms_pooled.mat');
    case {'contra', 'ipsi', 'contraipsi'};
        % load original crosscorrelated audiograms
        load('D:\Data\RueData\IBW\Audiograms\Corr_450ms\AllStats_450ms.mat');
    otherwise,
        error(['Invalid, illegal & utterly stupid keyword ''' kw '''.']);
end
% additional filtering
switch kw
    case 'contra',
        Stats = Stats(isRueContra(Stats));
    case 'ipsi',
        Stats = Stats(isRueIpsi(Stats));
end

% now we have a struct array Stats
ispair = isRuePair(Stats);
Pstats = Stats(ispair);
NPstats = Stats(~ispair);

Pstats = sortAccord(Pstats, [Pstats.mean]);
NPstats = sortAccord(NPstats, [NPstats.mean]);

Stats = [Pstats, NPstats];

% output 1: linear list, throw out autocorr
FN = ['LinList_' kw '.xls'];
L = Stats(~isRueDiag(Stats)); % no autocorrs
c1 = {L.cell_1}.';
c2 = {L.cell_2}.';
m = {L.mean}.';
s = {L.std}.';
MM = [c1 c2 m s];
xlswrite(FN, MM);

% output 2: matrix
FN1 = ['MatrixMeancorr_' kw '.xls'];
FN2 = ['MatrixStdcorr_' kw '.xls'];
FN3 = ['CellnameList_' kw '.xls'];
Name = unique({Stats.cell_1 Stats.cell_2});
Ncell = numel(Name);
Pnam = {Pstats.cell_1 Pstats.cell_2};
qp = ismember(Name, Pnam); % whether Name is one of a pair
ip = find(qp);
inp = find(~qp);
Pname = [sort(Name(ip)) sort(Name(inp))]; % sorted names
% visit all elements of Stats and put them in the appriate matrix element
Mmean = nan(Ncell);
Mstd = nan(Ncell);
for ii=1:numel(Stats),
    st = Stats(ii);
    i1 = strmatch(st.cell_1, Pname, 'exact');
    i2 = strmatch(st.cell_2, Pname, 'exact');
    Mmean(i1,i2) = st.mean;
    Mmean(i2,i1) = st.mean;
    Mstd(i1,i2) = st.std;
    Mstd(i2,i1) = st.std;
end
xlswrite(FN1, Mmean);
xlswrite(FN2, Mstd);
xlswrite(FN3, Pname');






