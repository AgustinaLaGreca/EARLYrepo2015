function s=nextto(G);
% GUIpiece/nextto - direction to place GUIpiece next to another GUIpiece.
%    nextto(G) returns the char string  'nextto Foo', where Foo is
%    name(G). This simplifies the input to GUIpiece/add, as in
%        Panel = PAdd(Panel, G1, nextto(G2)); 
%
%    See GUIpanel, GUIpiece/add, GUIpiece, below.

s = ['nextto ' name(G)];






