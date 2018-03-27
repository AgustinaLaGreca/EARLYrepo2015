function Y = getGUIdata(h, FN, Def);
% getGUIdata - get struct field of userdata property of a graphics object
%    getGUIdata(h, FN) returns the FN field of struct-valued 
%    'userdata' property of h. An error results if Fn is not a field
%    of the userdata.
%   
%    getGUIdata([h1 h2 ..], FN) returns multiple data in cell array.
%
%    getGUIdata(h) returns the entire struct in userdata.
%
%    getGUIdata(h, FN, DefaultValue) returns the FN field if it exists,
%    DefaultValue if not.
%
%    The names FN may also contain subfields, like 'Foo.fld', etc, in 
%    which case the 'fld' subfield of Foo is returned.
%
%    See also setGUIdata, rmGUIdata, findGUIdata, GUIgetter, GUIputter.


provideDefault = (nargin>2);

if length(h)>1, % recursive handling
    if nargin==1, Args = {};
    elseif nargin==2, Args = {FN};
    else, Args = {FN Def};
    end    
    for ii=1:length(h),
        Y{ii} = getGUIdata(h(ii), Args{:});
    end
    return
end

% handle FN subfields by recursion
if nargin>1,
    Mult_FN = Words2cell(FN,'.');
    if numel(Mult_FN)>1,
        Y = getGUIdata(h,Mult_FN{1}); % this delegates the error handling to the regular section of getGUIdata (below)
        % peel of the subfields, except last one
        for isub=2:numel(Mult_FN)-1,
            Y = Y.(Mult_FN{isub});
        end
        % handling of last subfield depends on availability of default value
        if provideDefault, 
            Y = getFieldOrDefault(Y,Mult_FN{end}, Def);
        else,
            Y = Y.(Mult_FN{end});
        end
        return;
    end
end
%----------single h, simple FN, from here------------

error(handleTest(h,'any', 'Graphics handle'));

ud = get(h,'userdata');
if isstruct(ud) & numel(ud)>1,
    error('Userdata contains struct array with multiple elements.')
end
if nargin<2,
    Y = ud;
    return;
end

FNfound = 0;
if isempty(ud),
elseif ~isfield(ud,FN),
else, 
    Y = ud.(FN);
    FNfound = 1;
end
    
if ~FNfound, 
    if provideDefault, Y = Def;
    else error(['Field ''' FN ''' not found in userdata.'])
    end
end





