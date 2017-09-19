function s=cornering(G);
% GUIpiece/cornering - direction to place a GUIpiece cornering another GUIpiece.
%    cornering(G) returns the char string  'cornering Foo', where Foo is
%    name(G). This simplifies the input to GUIpiece/add, as in
%        Panel = PAdd(Panel, G1, cornering(G2)); 
%
%    See GUIpanel, GUIpiece/add, GUIpiece/nextto, GUIpiece/cornering,
%    GUIpiece/ontopof.

s = ['cornering ' name(G)];






