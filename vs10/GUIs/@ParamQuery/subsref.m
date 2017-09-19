function Y = subsref(P, S);
% paramquery/SUBSREF - SUBSREF for paramquery objects
%
%   P.Foo is the the same as Foo(P) provided that Foo is a field of
%   struct(P) or a method for paramquery objects.
%
%   P(5) returns the 5th element of Paramquery array P.
%
%   P('Foo', 'Faa', ...) returns the subarray of array P containing the
%   elements named Foo, Faa, ... .
%
%   P('Foo*') returns the subarray of array P whose names start with Foo.
%
%   See also "methods paramquery".

if length(S)>1, % use recursion from the left
    y = subsref(P,S(1));
    Y = subsref(y,S(2:end));
    return;
end
%----single-level subsref from here (i.e., S is scalar) ---
switch S.type,
    case '()', % array element
        if iscellstr(S.subs), % look for named object
            nam = S.subs;
            Y = [];
            if isequal(1,numel(nam)) && ~isempty(nam{1}) && isequal('*', nam{1}(end)), % wildcards
                imatch = strmatch(nam{1}(1:end-1),{P.Name});
                Y = P(imatch);
            else, % exact matches
                for ii=1:numel(nam),
                    imatch = strmatch(nam{ii},{P.Name},'exact');
                    if ~isempty(imatch), Y = [Y P(imatch)]; end;
                end
            end
            if isequal({''}, nam) || isequal(cell(0,1), nam) || isequal(cell(1,0), nam), % return empty paramquery array
                Y = ParamQuery; Y = Y([]);
            elseif isempty(Y),
                error(['Array of ' class(P) ' object does not contain elements with the name(s) as specified in parentheses.'])
            end
        else,
            Y = builtin('subsref',P,S);
        end
    case '.',
        Field = S.subs; % field name
        [Fldname, Mess] = keywordMatch(Field, fieldnames(P), 'paramquery field name or method');
        if isempty(Mess), 
            Y = P.(Fldname);
        else, % try method
            [Method, Mess] = keywordMatch(Field, methods(P), 'paramquery field name or method');
            if isempty(Mess),
                Y = feval(Method, P);
            end
        end
        error(Mess);
    case '{}',
        error('paramquery objects do not allow subscripted reference using braces {}.');
end % switch/case
