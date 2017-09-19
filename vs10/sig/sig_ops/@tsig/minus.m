function P = minus(S,T);
% tsig/minus - subtract tsig objects
%   S-T is the channel-wise difference of tsig objects S and T. The channel 
%   counts of S and T must be compatible:  
%      [S.nchan T.nchan] == [M M] or [1 M] or [M 1], with M>=1.
%   The onsets of S and T are accounted for! The subtraction runs from 
%   min(S.onset,T.onset) to max(S.offset,T.offset) and S and T are aligned 
%   before being subtracted. In that part of the time domain that is inside S,
%   but outside T, S-T equals S by definition (and -T in the opposite case). 
%   This amounts to padding zeros in the uncovered parts of the time domains
%   of the two signals before subtracting them.
%
%   S-X, where S is a tsig object and X a scalar number, subtracts X from 
%   the samples of S. If X is a row vector, its elements are subtracted from
%   the corresponding channels of S. Note that subtracting a row vector from
%   a single-channel tsig has the side effect of creating a multi-channel
%   tsig.
%
%   X-S, where X is numerical and S a tsig object, equals -(S-X) by
%   definition. The same remarks apply as for S-X above.
%
%   See tsig/plus for the caveats involved when using multi-term 
%   expressions (without parenthesis) containing both tsig objects and 
%   numbers.
%
%   See also tsig/plus, tsig/uminus, sig_ops/@tsig.

% work is delegated to private function binop
[P, Mess] = binop(S,T,@minus,0,0);
% binop returns keyword-style messages and does not crash. Throw tailored error
switch Mess,
    case 'void', error('Void tsig objects may not be subtracted.');
    case 'numnonrowvector', error('Numerical values subtracted from tsig objects must be scalars or row vectors.');
    case 'nchan', error('Incompatible channel counts between sig objects or between tsig and numerical value.')
    case 'nontsig', error('Tsig objects can only be subtracted from numerical values or other tsig objects.');
    case 'fsam', error('Tsig object having different sample frequencies cannot be be subtracted.');
end

