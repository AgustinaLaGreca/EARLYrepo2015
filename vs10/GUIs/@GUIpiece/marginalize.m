function P=marginalize(P,XY);
% GUIpiece/marginalize - provide inner margins to GUIpiece
%   P=marginalize(P,[DX DY]) provides a right margin of DX pixels, and a
%   lower margin of DY pixels to GUIpiece P. This is realized by temporarily 
%   adding invisible children to P. Provide margins only after all children 
%   of P have been defined.
%
%    See also GUIpiece, GUIpiece/add, GUIpiece/draw.

V1 = struct('Name','q_void____1','Extent',[0.01 0.01]);
V2 = struct('Name','q_void____1','Extent',[0.01 0.01]);
P=add(P,V1,'nextto rightmost',XY);
P=add(P,V2,'below lowest',XY);
% now remove them - adding them has strechtched P's Extent, so the work is done
P.ChildArrangement = P.ChildArrangement(1:end-2);
P.Children = P.Children(1:end-2);



