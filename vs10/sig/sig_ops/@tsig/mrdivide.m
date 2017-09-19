function P = mrdivide(S,T);
% tsig/mrdivide - divide tsig objects
%   S/T or S./T is the channel-wise quotient of tsig objects S and T. 
%   The channel counts of S and T must be compatible:  
%      [S.nchan T.nchan] == [M M] or [1 M] or [M 1], with M>=1.
%   The onsets of S and T are accounted for! The division runs from 
%   min(S.onset,T.onset) to max(S.offset,T.offset) and S and T are aligned 
%   before being divided. In that part of the time domain which is outside S,
%   but inside T, S/T is set to 0 by definition. In the region outside T 
%   but inside S, S/T is set to S. This amounts to padding zeros and ones
%   in the out-of-time-domain values of S and T, respectively.
%
%   For tsig objects, there is no distinction between S./T and S/T.
%
%   S/X, where S is a tsig object and X a scalar number, divides the 
%   samples of S by X. If X is a row vector, its elements divide 
%   the corresponding channels of S. Note that dividing a single-channel 
%   tsig by a row vector has the side effect of creating a multi-channel
%   tsig.
%
%   X/S, where X is numerical and S a tsig object, equals X times the
%   reciprocal of S, 1/S. The same remarks on row vectors and channels
%   apply as for S/X above.
%
%   See tsig/plus for the caveats involved when using multi-term 
%   expressions (without parenthesis) consisting of both tsig objects and 
%   numbers.
%
%   See also tsig/plus, tsig/modulate, sig_ops/@tsig.

% work is delegated to private function binop
[P, Mess] = binop(S,T,@rdivide,0,1);
% binop returns keyword-style messages and does not crash. Throw tailored error
switch Mess,
    case 'void', error('Void tsig objects may not be divided.');
    case 'numnonrowvector', error('Numerical values divided by or dividing tsig objects must be scalars or row vectors.');
    case 'nchan', error('Incompatible channel counts between sig objects or between tsig and numerical value.')
    case 'nontsig', error('Tsig objects can only be divided by or divide numerical values or other tsig objects.');
    case 'fsam', error('Tsig object having different sample frequencies cannot be be divided.');
end



