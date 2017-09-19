function ip = isRuePair(M,M2)
% isRuePair - tRue for Rue pairs of the Stats output of RueCorrAllPairs

if nargin==2,
    n1 = M;
    n2 = M2;
elseif numel(M)>1,
    for ii=1:numel(M),
        ip(ii) = isRuePair(M(ii));
    end
    return;
else,
    n1 = M.cell_1;
    n2 = M.cell_2;
end



ip = ~isequal(n1, n2) && (...
    isequal(strrep(n1,'A', 'B'), n2) ...
    || isequal(strrep(n1,'B', 'A'), n2)   );




