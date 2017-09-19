function Y = cumsum(T, dim);
% tsig/cumsum - CUMSUM for tsig objects.
%    C=cumsum(S) for tsig object S is a tsig whose channel C(k) contains
%    the cumulative sum of the kth channel of S, S(k).
%  
%    cumsum(S,1) is the same as cumsum(S). By convention, the "first 
%    dimension" of tsig objects is the time dimension. This is consistent 
%    with tsig/vertcat appending tsig objects in time. To perform a proper
%    "time integral" of S, see tsig/trapz.
%
%    C=cumsum(S,2) returns a tsig object whose kth channel C(k) equals
%    S(1)+S(2)..S(k). Importantly, the channels are time-aligned
%    before being multiplied. Thus C(k) runs from min(T(1:k).onset) to
%    max(T(1:k).offset). 
%
%    To compute the indefinite time integral of T, see tsig/cumtrapz.
% 
%    See also tsig/trapz, tsig/dt, tsig/sum, tsig/plus.

if nargin<2, dim=1; end

if ~isequal(1,dim) && ~isequal(2,dim),
    error('Dim argument must be either 1 or 2 for tsig/cumsum.');
end
if isvoid(T),
    error('Invalid input: void tsig object.'); 
end

if dim==1, % cumsum of each channel
    Y = wavefun(@cumsum, T, 'tsig');
elseif dim==2, % mean waveform across channels: delegate to xchanop
    if nchan(T)==0,
        Y = T; return
    end
    Y = channel(T,1); % means Y=T(1), but cannot use such subsasgn syntax within tsig method
    for ii=2:nchan(T),
        Y = channel(Y, ii, channel(Y,ii-1) + channel(T,ii)); % means Y(ii)=Y(ii-1)+T(ii), but cannot use such subsasgn syntax within tsig method
    end
end



