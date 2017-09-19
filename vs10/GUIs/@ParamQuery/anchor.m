function A = anchor(Q);
% ParamQuery/anchor - anchoring point for horizontal alignment
%   anchor(Q) returns the width in pixels of the prompt of Q. This
%   coincides with the left edge of the edit part of Q, which is used for
%   horiontal algnment of paramqueries.
%
%   See also ParamQuery/add, ParamQuery/draw.

[dum S]=width(Q);
A = S.Prompt;




