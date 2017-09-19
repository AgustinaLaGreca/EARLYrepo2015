function P = plus(S,T);
% tsig/plus - add tsig objects
%   S+T is the channel-wise sum of tsig objects S and T. The channel counts
%   of S and T must be compatible:  
%      [S.nchan T.nchan] == [M M] or [1 M] or [M 1], with M>=1.
%   The onsets of S and T are accounted for! The sum runs from 
%   min(S.onset,T.onset) to max(S.offset,T.offset) and S and T are aligned 
%   before being added. In that part of the time domain that is inside S,
%   but outside T, S+T equals S by definition (and vice versa). This 
%   amounts to padding zeros in the uncovered parts of the time domains
%   of the two signals before adding them.
%
%   S+X==X+S, where S is a tsig object and X a scalar number, adds X to all
%   samples of S. If X is a row vector, its elements are added to the
%   corresponding channels of S. Note that adding a row vector to a
%   single-channel tsig has the side effect of creating a multi-channel
%   tsig.
%
%   Note that the alignment of signals before adding them, combined with
%   the option to add numbers to tsig objects, can turn the addition
%   into a non-associative operation. If S and T are tsig objects having
%   different onstes and/or offsets, and X is a number, then
%   (S+T)+X is not equal to S+(T+X). In the former case +X is applied to
%   the total time domain of (S+T), while in the latter case, +X is only
%   applied to the domain of T. This makes S+T+X ambiguous. Matlab R2007a,
%   interpretes S+T+X as (S+T)+X, but this might change. It is wise to
%   disambiguate multi-term sums by using ().
%
%   See also tsig/minus, tsig/uplus, sig_ops/@tsig.

% work is delegated to private function binop
[P, Mess] = binop(S,T,@plus,0,0);
% binop returns keyword-style messages and does not crash. Throw tailored error
switch Mess,
    case 'void', error('Void tsig objects may not be added.');
    case 'numnonrowvector', error('Numerical values added to tsig objects must be scalars or row vectors.');
    case 'nchan', error('Incompatible channel counts between sig objects or between tsig and numerical value.')
    case 'nontsig', error('Tsig objects can only be added to numerical values or other tsig objects.');
    case 'fsam', error('Tsig object having different sample frequencies cannot be be added.');
end
