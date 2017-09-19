function S = RueSelectSignCells(Twin, Selection);
% RueSelectSignCells - select cell pairs based on significance
if nargin<1, Twin = 450; end
if nargin<2, Selection = ''; end

% load diagonal cells & apply selection
DD = LookAtDiag(Twin);
if ~isempty(Selection),
    qok = [DD.(Selection)]==1;
    DD = DD(qok);
end
iwrong = strmatch('p090427Bi', {DD.cell_1}); 
[DD(iwrong).cell_1] = deal('p90427Bi');
[DD(iwrong).cell_2] = deal('p90427Bi');
OKnames = {DD.cell_1};

% read all pairs
FN = ['AllStats_' num2str(Twin) 'ms.mat'];
FFN = fullfile(['D:\Data\RueData\IBW\Audiograms\Corr_' num2str(Twin)  'ms'], FN);
qq = load(FFN);
CC = qq.Stats;
nam1 = {CC.cell_1}; nam2 = {CC.cell_2};
% fix incorrect naming
iwrong = strmatch('p090427Bi', nam1); [CC(iwrong).cell_1] = deal('p90427Bi');
iwrong = strmatch('p090427Bi', nam2); [CC(iwrong).cell_2] = deal('p90427Bi');
nam1 = {CC.cell_1}; nam2 = {CC.cell_2};

% both cells of a pair should belong to signifant group; filter for that
qok = ismember(nam1, OKnames) & ismember(nam2, OKnames);
CC = CC(qok);
nam1 = {CC.cell_1}; nam2 = {CC.cell_2};


allNames = unique(nam1);
% re-sort names: pairs first, then the rest
ipnames = strmatch('p', allNames); ienames = strmatch('e', allNames);
Npaired = numel(ipnames);
allNames = [sort(allNames(ipnames)) sort(allNames(ienames))];
Ncell = numel(unique(nam1));

% 

S = CollectInStruct(Twin, Selection, CC, allNames, Ncell, Npaired);




