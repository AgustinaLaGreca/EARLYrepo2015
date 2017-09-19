function imatch = find(Q, nam)
% PulldownMenu/find - find named CycleList element within array
%   find(P, 'Foo') returns the index or indices within PulldownMenu array
%   P of the element(s) named Foo. Matching is case sensitive and exact.
%   [] is returned if no match occurs.
%
%   See CycleList.

imatch = strmatch(nam, {Q.Name}, 'exact');





