function P=GUIpiece(Name, Content, Extent, OriginPos, LowRightMargins, Color);
% GUIpiece - construct GUIpiece object
%   P=GUIpiece(Name, Content, Extent, OriginPos, LowRightMargins, Color);

if nargin<3,
    Extent = [0 0];
end
if nargin<4,
    OriginPos = [0 0];
end
if nargin<5,
    LowRightMargins = [0 0];
end
if nargin<6,
    Color = [];
end

% specials cases: struct input or object input
if nargin==1 && isstruct(Name),
    P = Name;
    P = class(P, mfilename);
    return;
elseif nargin==1 && isa(Name, mfilename), % object in, object out
    P = Name;
    return;
elseif nargin==1, Content = [];
elseif nargin<1, % void object with correct fields
    [Name, Content, Children, ChildArrangement, Extent, Color] = deal([]);
    P = CollectInStruct(Name, Content, Children, ChildArrangement, LowRightMargins, Color, Extent);
    P = class(P,mfilename);
    return;
end

% generic case
if ~isvarname(Name),
    error('Name arg of GUIpiece constructor must be valid Matlab identifier.')
end
Children = {'Origin'};
ChildArrangement = struct('Name', 'Origin', 'RelPos',[],...
    'Neighbor',[],'Shift',[], 'Position', OriginPos, 'Extent', [0 0]);
P = CollectInStruct(Name, Content, Children, ChildArrangement, LowRightMargins, Color, Extent);
P = class(P,mfilename);




