function P=draw(h, P, XYpos);
% GUIpiece/draw - draw GUIpiece and its components
%   draw(hparent,P,XYpos) draws GUIpiece P in the graphics object with
%   handle hparent at position XYpos from the left-upper corner. Any 
%   children of P are drawn in a recursive way.
%
%    See also GUIpiece, GUIpiece/add.

if nargin<3, XYpos = [0 0]; end
if ~ishandle(h),
    error('Input argument h must be single graphics handle.');
end

if isequal('figure', get(h,'type')), % make sure P fits in figure
    set(h,'units','pixels');
    figpos = get(h,'position');
    lowerRight = P.Extent+XYpos;
    DX = max(0,lowerRight(1)-figpos(3));
    DY = max(0,lowerRight(2)-figpos(4));
    figpos = figpos+[0 -DY DX DY]; % adjust size of fig, make sure upper border stays put on screen
    set(h,'position',figpos);
end
if ~isempty(P.Content), % P may be an explicit uipanel
    ParentPos = get(h,'position');
    XYcontainer = localFigpos(XYpos, P.Extent, ParentPos(3:4)); % convert to Matlab fig coordinates w origin at lower-left corner
    h = uipanel('Parent', h, 'units', 'pixels', 'Position', [XYcontainer P.Extent], P.Content); % the explicit container will act as parent handle
    h = double(h);
    XYpos = [0 0]; % the shift of children is already realized by shifting the container
end
% draw children one by one
Ytop = P.Extent(2);
for ii=2:length(P.Children), % start at 2 - #1 is void 'Origin'
    c = P.Children{ii};
    ca = P.ChildArrangement(ii);
    if isstruct(c), % temp fake object for testing & debugging
        XY = localFigpos(XYpos+ca.Position, ca.Extent, P.Extent); % convert to Matlab fig coordinates w origin at lower-left corner
        uicontrol(h, 'units', 'pixels', 'position', [XY c.Extent], ...
            'style', 'text', 'string', c.Name, 'backgroundcolor', rand(1,3));
    else, % objects know how to draw themselves
        c=draw(h,c,ca.Position);
    end
end

if isGUI(h),
    GUIrepos(h);
    set(h,'visible','on');
end
%======================
function Pos = localFigpos(Pos, Ext, ContainerExtent);
% convert upper-left-based coordinates within container to lower-left-based ones
Pos = [Pos(1) ContainerExtent(2)-Pos(2)-Ext(2)];





