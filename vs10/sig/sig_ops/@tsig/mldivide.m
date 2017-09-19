function P = mldivide(S,T);
% tsig/ldivide - modulate tsig object by another tsig object
%   S\M, S.\M or modulate(S,M) modulates "carrier" S by multiplying S 
%   with "modulator" M. The result is the is the channel-wise product of 
%   tsig objects S and M. The channel counts of S and M must be compatible
%   (see tsig/plus).
%   The onsets of S and M are accounted for! The multiplication runs from 
%   min(S.onset,M.onset) to max(S.offset,M.offset) and S and M are aligned 
%   before being multiplied. In the region inside S but outside T, S\M is
%   equal to S; in the region outside S, S\M is zero. Note that these
%   conventions are subtly different from the product S*M. The conventions
%   for S\M can be described as "where M exists, M modulates S; where M 
%   does not exist, S is unaltered; where S does not exists, S\M is zero." 
%
%   For tsig objects, there is no distinction between S.\M and S\M.
%
%   S\X or X\S, where S is a tsig object and X a scalar number, multiplies 
%   the samples of S by X. This is comepletely equal to S*X X*S. See
%   tsig/times.
%
%   See tsig/plus for the caveats involved when using multi-term 
%   expressions (without parenthesis) consisting of both tsig objects and 
%   numbers.
%
%   See also tsig/times, tsig/plus, sig_ops/@tsig.

% work is delegated to private function binop
[P, Mess] = binop(S,T,@times,0,1); % the last 1 makes the difference re tsig/times
% binop returns keyword-style messages and does not crash. Throw tailored error
switch Mess,
    case 'void', error('Void tsig objects may not be multiplied.');
    case 'numnonrowvector', error('Numerical values modulating tsig objects must be scalars or row vectors.');
    case 'nchan', error('Incompatible channel counts between sig objects or between tsig and numerical value.')
    case 'nontsig', error('Tsig objects can only be be multiplied by numerical values or other tsig objects.');
    case 'fsam', error('Tsig object having different sample frequencies cannot be be multiplied.');
end



