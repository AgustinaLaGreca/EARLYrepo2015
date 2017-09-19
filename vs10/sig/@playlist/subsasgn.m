function P = subsasgn(P, S, B);
% playlist/SUBASGN - SUBSASGN for playlist objects
%
%   P.Foo=X  assigns value X to assignable property Foo of P.
%   Type P.help=0 to get a list of all fields that may be set by 
%   subscripted assigment
%
%   See also playlist, playlist/subsref, "methodshelp playlist".

if length(S)>1, % use recursion from the left.   P.Foo.Goo = X
    p = subsref(P,S(1));     %                   p = P.Foo
    p = subsasgn(p, S(2:end)); %                 p.Goo = X
    P = subsasgn(P,S(1),p); %                    P.Foo = p
    return;
end

%----single-level subsref from here (i.e., A is scalar) ---

switch S.type,
    case '()',
        error('Syntax P(..)=RHS is invalid for playlist objects.')
    case '.',
        AssignMethods = {'list' 'help'};
        Field = S.subs; % field name
        [Method, Mess] = keywordMatch(Field, AssignMethods, 'assigned playlist field name');
        error(Mess);
        if isequal('help', Method),
            methodsHelp(P, AssignMethods);
        else,
            P = feval(Method, P, B);
        end
    case '{}',
        error('Playlist objects do not allow subscripted reference using braces {}.');
end % switch/case


