function P = times(S,T);
% tsig/times - multiply tsig objects
%   S*T or S.*T is the channel-wise product of tsig objects S and T. 
%   The channel counts of S and T must be compatible:  
%      [S.nchan T.nchan] == [M M] or [1 M] or [M 1], with M>=1.
%   The onsets of S and T are accounted for! The multiplication runs from 
%   min(S.onset,T.onset) to max(S.offset,T.offset) and S and T are aligned 
%   before being multiplied. In the non-overlapping part of S and T, the 
%   product S*T is zero. This amounts to assuming that S and T are zero 
%   outside their respective domains. For multiplication using a different 
%   convention concerning the non-overlap region, see tsig/modulate.
%
%   For tsig objects, there is no distinction between S.*T and S*T.
%
%   S*X or X*S, where S is a tsig object and X a scalar number, multiplies 
%   the samples of S by X. If X is a row vector, its elements multiply 
%   the corresponding channels of S. Note that multiplying a single-channel 
%   tsig by a row vector has the side effect of creating a multi-channel
%   tsig.
%
%   See tsig/plus for the caveats involved when using multi-term 
%   expressions (without parenthesis) consisting of both tsig objects and 
%   numbers.
%
%   See also tsig/modulate, tsig/plus, sig_ops/@tsig.

% work is delegated to private function binop
[P, Mess] = binop(S,T,@times,0,0);
% binop returns keyword-style messages and does not crash. Throw tailored error
switch Mess,
    case 'void', error('Void tsig objects may not be multiplied.');
    case 'numnonrowvector', error('Numerical values multiplying tsig objects must be scalars or row vectors.');
    case 'nchan', error('Incompatible channel counts between sig objects or between tsig and numerical value.')
    case 'nontsig', error('Tsig objects can only be be multiplied by numerical values or other tsig objects.');
    case 'fsam', error('Tsig object having different sample frequencies cannot be be multiplied.');
end



