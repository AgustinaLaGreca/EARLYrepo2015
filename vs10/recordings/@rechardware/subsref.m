function Y = subsref(E, S);
% recHardware/SUBSREF - SUBSREF for recHardware objects
%
%   RH(5) is the 5th element of recHardware array E.
%
%   RH{..} is not a valid syntax.
%
%   RH.Foo returns properties of RH, i.e., fields of struct(RH). 
%
%   Staggered subsref usage as in RH.Foo(4).Faa is handled recursively.
%
%   See also dataset, "methods dataset".


if length(S)>1, % use recursion from the left
    y = subsref(E,S(1));
    Y = subsref(y,S(2:end));
    return;
end
%----single-level subsref from here (i.e., S is scalar) ---
switch S.type,
    case '()', % array element
        try, Y = builtin('subsref',E,S);
        catch, error(lasterr);
        end
    case '{}',
        error('Experiment objects do not allow subscripted reference using braces {}.');
    case '.',
        Field = S.subs; % field name
        Y = E.(Field);
end % switch/case


