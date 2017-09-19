function Y = haslist(P);
% playsig/haslist - true for playlist object that has its list specified.
%    haslist(P) or P.haslist returns true (1) if a list has been spefied 
%    for playlist T, false (0) otherwise.
%
%    See also playlist/list.

Y = ~isempty(P.iplay);





