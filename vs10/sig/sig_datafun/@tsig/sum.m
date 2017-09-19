function Y = sum(T, dim);
% tsig/sum - SUM for tsig objects.
%    sum(S) for tsig object S is a row vector M whise elements M(k) are 
%    the sum of all the samples in the respective channels of T, T(k).
%  
%    Sum(S,1) is the same as sum(S). Thus by convention, the "first 
%    dimension" of tsig objects is the time dimension. This is consistent 
%    with tsig/vertcat appending tsig objects in time.
%
%    T=sum(S,2) returns a tsig object containing the sum of the channels of
%    S: S = T(1)+T(2)+.. Thus by convention, the "second dimension" of tsig
%    objects is the channel number, which is consistent with tsig/horzcat
%    concatenating channels. Importantly, the channels are time-aligned
%    before being summed. The summed tsig T runs from min([S.onset]) to 
%    max(S.offset), and for channels that do not cover this full range, 
%    zeros are padded in the missing parts.
% 
%    See also tsig/max, tsig/mean, tsig/plus, tsig/cat, tsig/isvoid.

if nargin<2, dim=1; end

if ~isequal(1,dim) && ~isequal(2,dim),
    error('Dim argument must be either 1 or 2 for tsig/sum.');
end
if isvoid(T),
    error('Invalid input: void tsig object.'); 
end

if dim==1, % sum of each channel
    Y = wavefun(@sum, T, 'cat');
elseif dim==2, % mean waveform across channels: delegate to xchanop
    my_sum = @(x)sum(x,2); % pass 2nd arg = DIM
    [Y, Mess] = xchanop(T, my_sum, 0); % 0 is the out-of-time-range value; see help text
    error(Mess);
end



