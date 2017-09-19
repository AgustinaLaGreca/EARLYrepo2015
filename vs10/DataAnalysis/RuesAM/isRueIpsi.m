function ip = isRueIpsi(M)
% isRueIpsi - tRue for ipsilateral Rue pairs of the Stats output of RueCorrAllPairs


if numel(M)>1,
    for ii=1:numel(M),
        ip(ii) = isRueIpsi(M(ii));
    end
    return;
end

n1 = M.cell_1;
n2 = M.cell_2;

ip = isequal('i', lower(n1(end))) & isequal('i', lower(n2(end)));




