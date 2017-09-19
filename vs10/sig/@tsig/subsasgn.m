function P = subsasgn(P, S, B);
% tsig/SUBASGN - SUBSASGN for tsig objects
%
%   T(k)=S replaces channel k of tsig object T with single-channel tsig S. 
%   Note: if k>T.nchan, the number of channels is increased to k. Any
%   channels between T.nchan and k are set to [].
%
%   T.Foo=X  assigns value X to assignable property Foo of P.
%   Type P.help=0 to get a list of all fields that may be set by 
%   subscripted assigment
%
%   See also tsig, tsig/subsref, tsig/cut, "methodshelp tsig".

if length(S)>1, % use recursion from the left.   P.Foo.Goo = X
    p = subsref(P,S(1));     %                   p = P.Foo
    p = subsasgn(p, S(2:end),B); %               p.Goo = X
    P = subsasgn(P,S(1),p); %                    P.Foo = p
    return;
end

%----single-level subsref from here (i.e., S is scalar) ---

switch S.type,
    case '()', % set channels.  P(k)=B
        idx = S.subs;
        if length(idx)>1,
            error('Invalid multi-dimensional assignment of tsig T(k,..)=X.');
        end
        idx = idx{1};
        if isequal(':',idx), idx=1:nchan(P); end
        try,
            P = channel(P,idx,B);
        catch,
            error(lasterr);
        end
    case '.', % set property. P.Foo = B
        AssignMethods = {'onset' 'fsam' 'cut' 'help'};
        Field = S.subs; % field name
        [Method, Mess] = keywordMatch(Field, AssignMethods, 'assigned tsig field name');
        error(Mess);
        if isequal('help', Method),
            methodsHelp(P, AssignMethods);
        else,
            P = feval(Method, P, B);
        end
    case '{}',
        error('tsig objects do not allow subscripted reference using braces {}.');
end % switch/case


