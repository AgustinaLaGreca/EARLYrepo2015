function Y = subsref(P, S);
% tsig/SUBSREF - SUBSREF for tsig objects
%
%   T(K) or channel(T,K) is T restricted to the channels in index array K.
%
%   P.Foo is the the same as Foo(P) provided that Foo is a method for tsig
%   objects.
%
%   P.help produces a list of all fields that may be assessed by subscripted
%   reference of tsig objects. These "fields" are really methods, so
%   type help tsig/Foo for help on field Foo.
%
%   See also tsig, "methods tsig".


if length(S)>1, % use recursion from the left
    y = subsref(P,S(1));
    Y = subsref(y,S(2:end));
    return;
end

%----single-level subsref from here (i.e., S is scalar) ---

switch S.type,
    case '()',
        idx = S.subs;
        if length(idx)>2,
            error('Invalid indexing of playlist object. Too many indices');
        end
        K = idx{1}; % K selects channels (see help text)
        if length(idx)>1, % return samples in numerical matrix
            error('Multiple indexing T(i,j,..) is not allowed for tsig objects.')
        end
        Y = channel(P,K);
    case '.',
        Field = S.subs; % field name
        [Method, Mess] = keywordMatch(Field, methods(P), 'tsig field name or method');
        error(Mess);
        Y = feval(Method, P);
    case '{}',
        error('tsig objects do not allow subscripted reference using braces {}.');
end % switch/case


