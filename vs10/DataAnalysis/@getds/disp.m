function Str = disp(G);
% getds/displ - display Getds object
%   disp(G) displays Getds object G.
%   
%   Str = disp(G) returns string Str.

Str = sprintf('Getds object for experiment %s', G.XP);
if nargout<1,
    disp(Str);
end





