function ip = isRueContra(M)
% isRueContra - tRue for contralateral Rue pairs of the Stats output of RueCorrAllPairs


if numel(M)>1,
    for ii=1:numel(M),
        ip(ii) = isRueContra(M(ii));
    end
    return;
end

n1 = M.cell_1;
n2 = M.cell_2;

ip = isequal('c', lower(n1(end))) & isequal('c', lower(n2(end)));




