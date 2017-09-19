function T = sparsify(T, rdev);
% transfer/sparsify - economic approximation if Transfer object
%    sparsify(T, rtol) is a "decimated" version of T having a relative
%    deviation of at most rtol. See sparsify for details.
%
%    See Transfer, sparsify.

if numel(T)>1, % recursion
    for ii=1:numel(T),
        T(ii) = sparsify(T(ii), rdev);
    end
    return;
end

% single T from here

if isfilled(T),
    I = sparsify(T.Freq, T.Ztrf, rdev);
    T.Freq = T.Freq(I);
    T.Ztrf = T.Ztrf(I);
end











