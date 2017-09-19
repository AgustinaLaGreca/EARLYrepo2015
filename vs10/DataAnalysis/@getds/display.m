function Str = display(G);
% getds/display - display Getds object
%   display(G) displays Getds object G.
%   
%   Str = disp(G) returns string Str.

Str = sprintf('%s =\n     Getds object for experiment %s', inputname(1), G.XP);
if nargout<1,
    disp(Str);
end





