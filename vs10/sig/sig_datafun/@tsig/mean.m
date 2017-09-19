function Y = mean(T, dim,flag);
% tsig/mean - MEAN for tsig objects.
%    M=mean(S) for tsig object S is a row vector M whose elements M(k) are 
%    the mean values of the respective channels of T, T(k).
%  
%    Mean(S,1) is the same as mean(S). Thus by convention, the "first 
%    dimension" of tsig objects is the time dimension. This is consistent 
%    with tsig/vertcat appending tsig objects in time.
%
%    T=mean(S,2) returns a tsig object containing the sample-wise mean
%    across channels. Thus by convention, the "second dimension" of tsig
%    objects is the channel number, which is consistent with tsig/horzcat
%    concatenating channels. Importantly, the channels are time-aligned
%    before being averaged. The mean tsig T runs from min([S.onset]) to 
%    max(S.offset), and for channels that do not cover this full range, 
%    zeros are padded in the missing parts.
% 
%    M=mean(S,2,'nozeros') is the same as mean(S,2) except that no zeros are
%    padded. In a given time range, the mean will be evaluated across only
%    those channels that cover this range. To illustrate the difference, 
%    suppose that the only channel in S that covers time interval (t0,t1) is
%    the first channel, S(1). In this time inteval, mean(S,2,'nozeros') 
%    equals S(1), while mean(S,2) equals S(1)/S.nchan.
%
%    See also tsig/max, tsig/std, tsig/plus, tsig/exist, tsig/isvoid.

if nargin<2, dim=1; end
if nargin<3, flag=''; end

if ~isequal(1,dim) && ~isequal(2,dim),
    error('Dim argument must be either 1 or 2 for tsig/mean.');
end
if isvoid(T),
    error('Invalid input: void tsig object.'); 
end
if ~isempty(flag),
    [flag, Mess] = keywordMatch(flag, {'nozeros'}, 'flag');
    error(Mess);
end

if dim==1, % mean of each channel
    Y = wavefun(@median, T, 'cat');
elseif dim==2, % mean waveform across channels: delegate to xchanop
    % 'nozeros' flag determines the default value passed to xchanop
    if isequal('nozeros',flag),
        defS = 'existingchannelsonly'; % exclude non-existing channel intervals
    else,
        defS = 0; % substitute zeros where channels don't exist
    end
    my_mean = @(x)mean(x,2); % pass 2nd arg = DIM
    [Y, Mess] = xchanop(T, my_mean, defS);
    error(Mess);
end



