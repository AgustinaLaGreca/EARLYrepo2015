function [Y,Tm] = min(T, T2);
% tsig/min - MIN for tsig objects.
%    M=min(S) for tsig object S is a row vector M with M(k) the minimum
%    sample value in S(k), the kth channel S.
%
%    [M,Tm] = min(S) also returns returns the time values in a row vector
%    Tm, with Tm(k) the time at which the minimum occurred in channel k.
%    If the miniumum value occurs more than once within a channel, the
%    first occurrence is returned.
%
%    min(S,T), where both S and T are tsig objects, is a tsig object
%    containing the channelwise and samplewise minimum among S and T.
%    S and T must have compatible channel; sizes (see tsig/plus).
%    S and T are time-aligned before being compared! In the time range
%    covered by S, but not by T, min(S,T) equals S, and vice versa. This
%    amounts to setting the value of the signals to inf outside their own
%    time range.
%
%    min(S,X) or min(X,S), with X a numeric value, is a tsig containing
%    the minimum of the sampes of S and the value X. For a row vector X,
%    the minimum is evaluated per channel. Note that specifying a
%    single-channel S and a row vector X turns the result into a
%    length(X)-channel tsig.
%
%    Any NaN's are ignored. It is an error to pass void tsig objects to
%    tsig/min.
%
%    See also min, tsig/min, tsig/datafun, tsig/isvoid.

if nargin==1,
    if isvoid(T), error('Invalid input: void tsig object.'); end
    [Y,I] = wavefun(@min, T, 'cat');
    Tm = onset(T)+DT(T)*(I-1); % convert sample indices to time
elseif nargin==2,
    if nargout>1, error('Max with two inputs does not return a second output.'); end
    [Y, Mess] = binop(T,T2,@min,inf,inf);
    switch Mess,
        case 'void', error('Invalid input: void tsig object.');
        case 'nchan', error('Incompatible channel counts.');
        case {'numnonrowvector', 'nontsig'}, error('Any non-tsig input argument must be numerical row vector.');
        case 'fsam', error('Unequal sample frequencies among tsig objects.');
    end
end



