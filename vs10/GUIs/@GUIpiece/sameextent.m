function P=sameExtent(P,Q,Dim);
% GUIpiece/sameExtent - match extents of GUIpiece
%   P=sameExtent(P,Q,'X') adjusts the width of P to match that of Q.
%   P=sameExtent(P,Q,'Y') adjusts the height of P to match that of Q.
%   P=sameExtent(P,Q,'XY') adjusts the width & height of P to match those of Q.
%
%    See also GUIpiece, GUIpanel.

[Dim, Mess] = keywordMatch(Dim, {'X' 'Y' 'XY'}, 'Dim input argument ');
error(Mess);

switch Dim,
    case 'X',
        P.Extent(1) = Q.Extent(1);
    case 'Y',
        P.Extent(2) = Q.Extent(2);
    case 'XY',
        P.Extent = Q.Extent;
end

