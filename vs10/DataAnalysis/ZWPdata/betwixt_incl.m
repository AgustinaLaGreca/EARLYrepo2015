function ib = betwixt_incl(X, A, B);
% BETWIXT - test if value is between two given values are equal to them
%    BETWIXT(X,A,B) or BETWIXT(X, [A B]) returns true for those elements
%    X(k) of X  for which A<=X(k)<=B.
%
%    For a combination of array and matrix inputs for X,A,B, these three
%    input arguments are "SameSized" prior to testing the inequality. In
%    such "mixed cases", it is best to use the Betwixt_incl(X,A,B) syntax,
%    as opposed to the Betwixt_incl(X,[A B]) syntax, which can be confusing.
%
%    See also betwixt.


error(nargchk(2,3,nargin));
if nargin<3,
   B = A(:,end);
   A = A(:,1);
end
if isempty(X), 
    ib = false(size(X)); 
else,
    [X,A,B] = SameSize(X,A,B);
    ib = (X>=A) & (X<=B);
end






