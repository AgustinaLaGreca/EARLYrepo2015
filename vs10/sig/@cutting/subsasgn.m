function U=subsasgn(U,S,B);
% subsasgn - subsasgn for cutting objects.
%    U(par1,par2,, .., B) performs a cut-and-replace on the object X 
%    contained in U. That is,
%
%        U(par1,par2, .. ,B) = cut(U.X,par1,par2, .., cutting(B)), 
%
%    which is X with the segment indicated by par1,par2 replaced by B.
%    The details of the cut-and-replace operation depend on the cut method
%    for class(U.X).
%
%    See also tsig/cut.

if length(S)>1, % use recursion from the left
    error(['Invalid staggering of subscripted assignment.' char(10) ...
        'Cut qualifiers can be followed by at most one () reference, as in A.cut(..) = X.'])
    return;
end

%----single-level subsref from here (i.e., S is scalar) ---
if isequal('.', S.type),
    error('Invalid cutting syntax X.cut.foo. Use X.cut(...) = ...')
end
if isequal('{}', S.type),
    error('Invalid cutting syntax X.cut{}. Use X.cut(...) = ...')
end
% syntax was X.cut(par1,..)
idx = S.subs; % index in cell array
U = cut(U.X, idx{:}, cutting(B)); % the class of the last object tells class(U.X)/cut that we're dealing with a cut-and-replace call



