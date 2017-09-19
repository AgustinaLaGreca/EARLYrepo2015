function Y = prod(T, dim);
% tsig/prod - PROD for tsig objects.
%    prod(S) for tsig object S is a row vector M whise elements M(k) are 
%    the products of all the samples in the respective channels of T, T(k).
%  
%    prod(S,1) is the same as prod(S). Thus by convention, the "first 
%    dimension" of tsig objects is the time dimension. This is consistent 
%    with tsig/vertcat appending tsig objects in time.
%
%    T=prod(S,2) returns a tsig object containing the product of the channels of
%    S: S = T(1)+T(2)+.. Thus by convention, the "second dimension" of tsig
%    objects is the channel number, which is consistent with tsig/horzcat
%    concatenating channels. Importantly, the channels are time-aligned
%    before being multiplied. The multiplied tsig T runs from min([S.onset]) to 
%    max(S.offset), and for channels that do not cover this full range, 
%    ones are padded in the missing parts.
% 
%    See also tsig/max, tsig/mean, tsig/plus, tsig/cat, tsig/isvoid.

if nargin<2, dim=1; end

if ~isequal(1,dim) && ~isequal(2,dim),
    error('Dim argument must be either 1 or 2 for tsig/prod.');
end
if isvoid(T),
    error('Invalid input: void tsig object.'); 
end

if dim==1, % prod of each channel
    Y = wavefun(@prod, T, 'cat');
elseif dim==2, % mean waveform across channels: delegate to xchanop
    my_prod = @(x)prod(x,2); % pass 2nd arg = DIM
    [Y, Mess] = xchanop(T, my_prod, 1); % 1 is the out-of-time-range value; see help text
    error(Mess);
end



