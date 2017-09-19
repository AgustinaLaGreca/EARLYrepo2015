function [Y, Tm] = max(T, T2);
% tsig/max - MAX for tsig objects.
%    M=max(S) for tsig object S is a row vector M with M(k) the maximum
%    sample value in S(k), the kth channel S.
%  
%    [M,Tm] = max(S) also returns returns the time values in a row vector
%    Tm, with Tm(k) the time at which the maximum occurred in channel k. 
%    If the maxiumum value occurs more than once within a channel, the 
%    first occurrence is returned.
% 
%    MAX(S,T), where both S and T are tsig objects, is a tsig object
%    containing the channelwise and samplewise maximum among S and T.
%    S and T must have compatible channel; sizes (see tsig/plus).
%    S and T are time-aligned before being compared! In the time range
%    covered by S, but not by T, max(S,T) equals S, and vice versa. This
%    amounts to setting the value of the signals to -inf outside their own
%    time range.
%
%    MAX(S,X) or MAX(X,S), with X a numeric value, is a tsig containing
%    the maximum of the sampes of S and the value X. For a row vector X,
%    the maximum is evaluated per channel. Note that specifying a
%    single-channel S and a row vector X turns the result into a
%    length(X)-channel tsig.
%
%    Any NaN's are ignored. It is an error to pass void tsig objects to
%    tsig/max.
%
%    See also MAX, tsig/min, tsig/datafun, tsig/isvoid.

if nargin==1,
    if isvoid(T), error('Invalid input: void tsig object.'); end
    [Y, I] = wavefun(@max, T, 'cat');
    Tm = onset(T)+DT(T)*(I-1); % convert sample indices to time
elseif nargin==2,
    if nargout>1, error('Max with two inputs does not return a second output.'); end
    [Y, Mess] = binop(T,T2,@max,-inf,-inf);
    switch Mess,
        case 'void', error('Invalid input: void tsig object.');
        case 'nchan', error('Incompatible channel counts.');
        case {'numnonrowvector', 'nontsig'}, error('Any non-tsig input argument must be numerical row vector.');
        case 'fsam', error('Unequal sample frequencies among tsig objects.');
    end
end



