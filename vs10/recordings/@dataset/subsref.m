function Y = subsref(P, S);
% dataset/SUBSREF - SUBSREF for dataset objects
%
%   DS(5) is the 5th element of dataset array DS.
%
%   DS.Foo returns the Foo field of struct(DS). If no field Foo exists,
%   DS.Foo is the the same as Foo(P) provided that Foo is a method for 
%   dataset objects. If Foo is neither a field nor a method of DS, DS.Foo
%   equals DS.Data.Foo, if Foo is a field of DS.Data. The last possibility 
%   to make sense of DS.Foo is DS.Stim.Foo. In short (denoting struct(DS)
%   as ds), DS.Foo is the first possible of the following:
%        ds.Foo       (field)
%        Foo(DS)      (method)
%        ds.Data.Foo  (field of ds.Data)
%        ds.Stim.Foo  (field of ds.Stim)
%   Staggered subseref usage as in DS.Foo(4).Faa is handled recursively.
%
%   DS.help displays this help text.
%
%   See also dataset, "methods dataset".


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
        DescrStr = 'dataset field name, method, or subfield';
        if isequal('help', Field), help('dataset/subsref'); return; end
        % try fieldname
        [Fldname, Mess] = keywordMatch(Field, fieldnames(P), DescrStr);
        if isempty(Mess), 
            Y = P.(Fldname);
        else, % try method
            [Method, Mess] = keywordMatch(Field, methods(P, '-full'), DescrStr);
            if isempty(Mess),
                Y = feval(Method, P);
            else, % try P.Data fieldname
                [Fldname, Mess] = keywordMatch(Field, fieldnames(P.Data), DescrStr);
                if isempty(Mess),
                    Y = P.Data.(Fldname);
                else, % try P.Stim fieldname
                    [Fldname, Mess] = keywordMatch(Field, fieldnames(P.Stim), DescrStr);
                    if isempty(Mess),
                        Y = P.Stim.(Fldname);
                    else,
                    end
                end
            end
        end
        error(Mess);
    case '{}',
        error('dataset objects do not allow subscripted reference using braces {}.');
end % switch/case
