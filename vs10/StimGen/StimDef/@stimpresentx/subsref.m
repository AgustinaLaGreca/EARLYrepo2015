function Y = subsref(P, S);
% StimPresent/SUBSREF - SUBSREF for StimPresent objects
%
%   P(5) is the 5th element of StimPresent array P.
%
%   P.Foo is the the same as Foo(P) provided that Foo is a method for 
%   StimPresent objects.
%
%   See also StimPresent, "methods StimPresent".

if length(S)>1, % use recursion from the left
    y = subsref(P,S(1));
    Y = subsref(y,S(2:end));
    return;
end
%----single-level subsref from here (i.e., S is scalar) ---
switch S.type,
    case '()', % array element
        Y = builtin('subsref',P,S);
    case '.',
        Field = S.subs; % field name
        if isfield(struct(P), Field),
            Y = P.(Field);
        else, %try abbreviated or none-case-matching field, or method
            [Fldname, Mess] = keywordMatch(Field, fieldnames(P), 'StimPresent field name or method');
            if isempty(Mess),
                Y = P.(Fldname);
            else, % try methods
                [Method, Mess] = keywordMatch(Field, methods(P), 'StimPresent field name or method');
                if isempty(Mess),
                    Y = feval(Method, P);
                end
            end
            error(Mess);
        end
    case '{}',
        error('StimPresent objects do not allow subscripted reference using braces {}.');
end % switch/case
