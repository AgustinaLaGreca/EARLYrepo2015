function P=add(P,p,PosInfo,Shift);
% GUIpiece/add - add component to GUIpiece
%   P=add(P,p,PosInfo,Shift) adds a component p to GUIpiece P.
%   p can be a GUIpiece or elementary component; p must have a Name and 
%   Extent. PosInfo is a char string like 'below Ref', where Ref is the 
%   name of a component previously added to P, or the string 'Origin', 
%   which denotes the upper-left corner of P, or one of
%     'last': last child added to P
%     'rightmost': rightmost positioned child 
%     'lowest': lowest positioned child 
%   Default is 'last'
%
%   The first word of PosInfo indicating the position re the Ref is one of
%     'below': below Ref, same left corner
%     'alignedwith': below Ref with their anchors alined. Anchor(p) and
%                Anchor(Ref) must exist (e.g. ParamQueries allow aligment).
%     'nextto': to the right of Ref, same vertical position
%     'cornering': upper-left corner of p touches lower-right corner of Ref
%     'ontopof': upper-left corners of Ref and p coincide.
%
%   Shift is an optional 2-component array [Xshift Yshift] containing an 
%   additional position shift in pixels. (Yshift is downward!)
%
%   Example
%     P=GUIpiece('FSmenu') % the complete GUI
%     P=add(P,FreqPanel); % freq params in upper-left corner
%     P=add(P, SPL, 'below FreqPanel', [10 0]); % below freq panel, 10 pixels to the right
%
%    See also GUIpiece, GUIpiece/draw.

if nargin<3, PosInfo = 'below'; end
if nargin<4, Shift = [0 0]; end

if nargout<1, 
    error('GUIpiece/add must be called using an output argument');
end

firstAdded = isequal(1,numel(P.Children)); % the origin is only one -> no true components added yet
if firstAdded && isa(p, 'ParamQuery'), % add separator indicating P's name to p. Used by GUIval ...
    p = separator(p, P.Name);  % ... for prettyfying its output struct
end

% figure out what the position of p is
PosInfo=Words2cell(PosInfo);
if length(PosInfo)<2, % previous one added (or Origin in case of first one)
    PosInfo{2} = 'last';
end
% first word is qualifier like below, etc
[Direction, Mess] = keywordMatch(PosInfo{1}, {'below' 'alignedwith' 'nextto' 'cornering' 'ontopof'},'direction specifier');
error(Mess);
% second reference child, either by name or location
Reference = localFindRef(PosInfo{2}, P.ChildArrangement); % decode locations like 'last', 'rightmost', etc
[Reference, Mess]=keywordMatch(Reference,{P.ChildArrangement.Name},'Component of GUIpanel');
error(Mess);

iref = strmatch(Reference,{P.ChildArrangement.Name}); % index of reference child
refPos = P.ChildArrangement(iref).Position; % position of reference child
refExt = P.ChildArrangement(iref).Extent; % extent of reference child
refch = P.Children{iref}; % Referent child
switch Direction,
    case 'below', newpos = [refPos(1) refPos(2)+refExt(2)] + Shift;
    case 'alignedwith', newpos = [refPos(1)+anchor(refch)-anchor(p) refPos(2)+refExt(2)] + Shift;
    case 'nextto', newpos = [refPos(1)+refExt(1) refPos(2)] + Shift;
    case 'cornering', newpos = [refPos(1)+refExt(1) refPos(2)+refExt(2)] + Shift;
    case 'ontopof', newpos = [refPos(1) refPos(2)] + Shift;
end
% add child & info to P
newChArr = struct('Name', p.Name, 'RelPos', Direction, 'Neighbor', Reference, ...
    'Shift', Shift, 'Position', newpos, 'Extent', p.Extent);
P.ChildArrangement(end+1) = newChArr;
P.Children{end+1} = p;
% update total extent of P
BottomRight = max(cat(1,P.ChildArrangement.Position) + cat(1,P.ChildArrangement.Extent));
XYorigin = P.ChildArrangement(1).Position; % use position of "origin" also to define lower-right margins
P.Extent = max(P.Extent, BottomRight+P.LowRightMargins);

%=================================================================
function Referent = localFindRef(Qualifier, CA)
Referent = Qualifier; % default: name of (hopefully) existing child
[Qualifier, Mess] = keywordMatch(Qualifier, {'last', 'rightmost', 'lowest'});
if ~isempty(Mess), % not a valid qualifier. Let calling block decide whether valid reference
    return;
end
pos = cat(1,CA.Position);
ext = cat(1,CA.Extent);
ilast = length(CA);
[~, ilowest] = max(pos(:,2)+ext(:,2));
[~, irightmost] = max(pos(:,1)+ext(:,1));
switch Qualifier,
    case 'last', iref = ilast;
    case 'rightmost', iref = irightmost;
    case 'lowest', iref = ilowest;
end
Referent = CA(iref).Name;


