function Y=subsref(U,S);
% subsref - subsref for cutting objects.
%    U(par1,...) performs a cut on the object X contained in U. That is,
%    U(par1,...) = cut(X,par1,...);
%
%    See also tsig/cut.

if length(S)>1, % use recursion from the left
    y = subsref(U,S(1));
    Y = subsref(y,S(2:end));
    return;
end

%----single-level subsref from here (i.e., S is scalar) ---
if isequal('.', S.type),
    error('Invalid cutting syntax X.cut.foo. Use X.cut(...).')
end
if isequal('{}', S.type),
    error('Invalid cutting syntax X.cut{}. Use X.cut(...).')
end
% syntax was X.cut(par1,..)
idx = S.subs; % index in cell array
Y = cut(U.X, idx{:});
