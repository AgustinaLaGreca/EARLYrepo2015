function P = mpower(S,T);
% tsig/mpower - power raising of tsig objects
%   S^T or S.^T is the channel-wise realization of S to the power T.
%   The channel counts of S and T must be compatible:  
%      [S.nchan T.nchan] == [M M] or [1 M] or [M 1], with M>=1.
%   The onsets of S and T are accounted for! The power raising runs from 
%   min(S.onset,T.onset) to max(S.offset,T.offset) and S and T are aligned 
%   before being combined. In that part of the time domain which is outside S,
%   but inside T, S^T is set to 0 by definition. In the region outside T 
%   but inside S, S^T is set to S. This amounts to padding zeros and ones
%   in the out-of-time-domain values of S and T, respectively.
%
%   For tsig objects, there is no distinction between S.^T and S^T.
%
%   S^X, where S is a tsig object and X a scalar number, raises the 
%   samples of S to the power X. If X is a row vector, its elements are the
%   powers of corresponding channels of S. Note that power raising a 
%   single-channel tsig by a row vector has the side effect of creating a 
%   multi-channel tsig.
%
%   X^S, where X is numerical and S a tsig object, equals X to the power 
%   the samples of S. The same remarks on row vectors and channels
%   apply as for S^X above.
%
%   See tsig/plus for the caveats involved when using multi-term 
%   expressions (without parenthesis) consisting of both tsig objects and 
%   numbers.
%
%   See also tsig/plus, tsig/uminus, sig_ops/@tsig.

% work is delegated to private function binop
[P, Mess] = binop(S,T,@power,0,1);
% binop returns keyword-style messages and does not crash. Throw tailored error
switch Mess,
    case 'void', error('Void tsig objects may not be power raised.');
    case 'numnonrowvector', error('Numerical values involved in power raising of tsig objects must be scalars or row vectors.');
    case 'nchan', error('Incompatible channel counts between sig objects or between tsig and numerical value.')
    case 'nontsig', error('Power raising of Tsig objects can only involve numerical values or other tsig objects.');
    case 'fsam', error('Tsig object having different sample frequencies cannot be be combined by power raising.');
end



