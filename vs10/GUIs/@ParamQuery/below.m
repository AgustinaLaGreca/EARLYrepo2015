function s=below(G);
% GUIpiece/below - direction to place GUIpiece below another GUIpiece.
%    below(G) returns the char string  'below Foo', where Foo is
%    name(G). This simplifies the input to GUIpiece/add, as in
%        Panel = PAdd(Panel, G1, below(G2)); 
%
%    See GUIpanel, GUIpiece/add, GUIpiece, nextto, alignedwith.

s = ['below ' name(G)];






