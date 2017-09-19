function Y = cumprod(T, dim);
% tsig/cumprod - CUMPROD for tsig objects.
%    C=cumprod(S) for tsig object S is a tsig whose channel C(k) contains 
%    the cumulative product of the kth channel of S, S(k).
%  
%    cumprod(S,1) is the same as cumprod(S). By convention, the "first 
%    dimension" of tsig objects is the time dimension. This is consistent 
%    with tsig/vertcat appending tsig objects in time. 
%
%    C=cumprod(S,2) returns a tsig object whose kth channel C(k) equals
%    S(1)*S(2)*..S(k). Importantly, the channels are time-aligned
%    before being multiplied. Zeros are padded outside ranges of existence.
%    Thus C(k) runs from min(T(1:k).onset) to max(T(1:k).offset), but it
%    nonzero only in the range max(T(1:k).onset) to min(T(1:k).offset)
% 
%    See also tsig/cumsum, tsig/dt, tsig/prod.

if nargin<2, dim=1; end

if ~isequal(1,dim) && ~isequal(2,dim),
    error('Dim argument must be either 1 or 2 for tsig/cumprod.');
end
if isvoid(T),
    error('Invalid input: void tsig object.'); 
end

if dim==1, % cumprod of each channel
    Y = wavefun(@cumprod, T, 'tsig');
elseif dim==2, % mean waveform across channels: delegate to xchanop
    if nchan(T)==0,
        Y = T; return
    end
    Y = channel(T,1); % means Y=T(1), but cannot use such subsasgn syntax within tsig method
    for ii=2:nchan(T),
        Y = channel(Y, ii, channel(Y,ii-1) * channel(T,ii)); % means Y(ii)=Y(ii-1)*T(ii), but cannot use such subsasgn syntax within tsig method
    end
end



