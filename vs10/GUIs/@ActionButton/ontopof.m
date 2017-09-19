function s=ontopof(G);
% GUIpiece/ontopof - direction to place a GUIpiece on top of another GUIpiece.
%    ontopof(G) returns the char string  'ontopof Foo', where Foo is
%    name(G). This simplifies the input to GUIpiece/add, as in
%        Panel = PAdd(Panel, G1, ontopof(G2)); 
%
%    See GUIpanel, GUIpiece/add, GUIpiece/nextto, GUIpiece/ontopof,
%    GUIpiece/ontopof.

s = ['ontopof ' name(G)];






