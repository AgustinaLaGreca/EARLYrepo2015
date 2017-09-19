function s=alignedwith(G);
% GUIpiece/alignedwith - direction to align a GUIpiece with another GUIpiece.
%    alignedwith(G) returns the char string  'alignedwith Foo', where Foo is
%    name(G). This simplifies the input to GUIpiece/add, as in
%        Panel = PAdd(Panel, G1, alignedwith(G2)); 
%
%    See GUIpanel, GUIpiece/add, GUIpiece/nextto, GUIpiece/cornering,
%    GUIpiece/ontopof.

s = ['alignedwith ' name(G)];






