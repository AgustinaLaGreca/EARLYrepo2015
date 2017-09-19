function imatch = find(Q, nam)
% CycleList/find - find named CycleList element within array
%   find(CL, 'Foo') returns the index or indices within cycleList array CL 
%   of the element(s) named Foo. Matching is case sensitive and exact.
%   [] is returned if no match occurs.
%
%   See CycleList.

imatch = strmatch(nam, {Q.Name}, 'exact');





