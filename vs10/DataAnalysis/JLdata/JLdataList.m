function LL=JLdataList;
%   JLdataList - list all recordings of the cells in JLbeat list
%   L = JLdatasearch('Foo', iCell) returns list L of all recordings from the
%   cell # iCell of experiment Foo.
%
%   See also JLdataList
L = {};
ML = JLreadBeats('list');
ExpIDs = fieldnames(ML);
for ii=1:numel(ExpIDs),
    expid = ExpIDs{ii};
    for icell = ML.(expid)(:).',
        [dum cs] = JLdatasearch(expid, icell);
        L = [L; cs];
    end
end
% sort
Nds = size(L,1)
sortParam = zeros(1,Nds);
for ii=1:Nds,
    xp = L{ii,1};
    iexp = str2double(xp(3:7));
    rc = L{ii,2};
    idash = find(rc=='-');
    icell = str2double(rc(2:idash(1)-1));
    irec = str2double(rc(idash(1)+1:idash(2)-1));
    sortParam(ii) = 1e6*iexp+1e2*icell+irec;
end
[dum, isort] = sort(sortParam);
L = L(isort,:);

LL = '';
for icol = 1:size(L,2),
    LL = [LL, strvcat(L(:, icol))];
end



% {JLp.SGSRExpID};
% {JLp.SGSRidentifier};
%SGSRdataset(JLp(1).SGSRExpID, JLp(1).SGSRidentifier)









