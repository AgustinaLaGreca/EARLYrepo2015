function Y = median(T, dim,flag);
% tsig/median - MEDIAN for tsig objects.
%    M=median(S) for tsig object S is a row vector M whose elements M(k) are 
%    the median values of the respective channels of T, T(k).
%  
%    Median(S,1) is the same as mean(S). Thus by convention, the "first 
%    dimension" of tsig objects is the time dimension. This is consistent 
%    with tsig/vertcat appending tsig objects in time.
%
%    T=median(S,2) returns a tsig object containing the sample-wise median
%    across channels. Thus by convention, the "second dimension" of tsig
%    objects is the channel number, which is consistent with tsig/horzcat
%    concatenating channels. Importantly, the channels are time-aligned
%    before being processed. The median tsig T runs from min([S.onset]) to 
%    max(S.offset), and for channels that do not cover this full range, 
%    zeros are padded in the missing parts.
% 
%    M=median(S,2,'nozeros') is the same as median(S,2) except that no zeros 
%    are padded. In a each given time range, the evaluation of the median is 
%    restricted to only those channels that exist during this range. See tsig/mean 
%    for more explanation.
%
%    See also tsig/mean, tsig/exist, tsig/plus, tsig/cat, tsig/isvoid.

if nargin<2, dim=1; end
if nargin<3, flag=''; end

if ~isequal(1,dim) && ~isequal(2,dim),
    error('Dim argument must be either 1 or 2 for tsig/median.');
end
if isvoid(T),
    error('Invalid input: void tsig object.'); 
end
if ~ischar(flag),
    error('Third argument must be character string.');
end
if ~isempty(flag),
    [flag, Mess] = keywordMatch(flag, {'nozeros'}, 'flag');
    error(Mess);
end

if dim==1, % median of each channel
    Y = wavefun(@median, T, 'cat');
elseif dim==2, % median waveform across channels: delegate to xchanop
    % 'nozeros' flag determines the default value passed to xchanop
    if isequal('nozeros',flag),
        defS = 'existingchannelsonly'; % exclude non-existing channel intervals
    else,
        defS = 0; % substitute zeros where channels don't exist
    end
    mymedian = @(x)median(x,2); % pass 2nd arg = DIM
    [Y, Mess] = xchanop(T, mymedian, defS);
    error(Mess);
end



