function P=remove(P, Name);
% GUIpiece/remove - remove component from GUIpiece
%   P = remove(P,'Foo') removes component named Foo from GUIpiece P.
%   This is useful for omitting certain paramQuery from a GUI without
%   having to rearrange the other queries.
%
%   P = remove(P, Q) where Q is an object (e.g. a ParamQuery) is the same
%   as remove(P, name(Q)).
%
%   See also GUIpiece/add.

if isempty(Name),
    return;
elseif ischar(Name),
    Name = cellify(Name);
elseif ~iscellstr(Name),
    p = Name;
    Name = {};
    for ii=1:numel(p),
        Name{ii} = name(p(ii));
    end
end

% now Name is a cellstring
for ii=1:numel(Name),
    nam = Name{ii};
    ihit = strmatch(nam, {P.ChildArrangement.Name}, 'exact');
    if isempty(ihit),
        error(['No query named ''' nam ''' in Guipiece.'])
    end
    P.Children(ihit) = [];
    P.ChildArrangement(ihit) = [];
end










