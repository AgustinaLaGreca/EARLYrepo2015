function ip = isRueDiag(M)
% isRueDiag - tRue for comparisons of a cell with itselfRue in the Stats output of RueCorrAllPairs


if numel(M)>1,
    for ii=1:numel(M),
        ip(ii) = isRueDiag(M(ii));
    end
    return;
end

n1 = M.cell_1;
n2 = M.cell_2;

ip = isequal(n1, n2);




