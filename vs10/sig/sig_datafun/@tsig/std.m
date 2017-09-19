function Y = std(T, dim,flag);
% tsig/std - STD for tsig objects.
%    M=std(S) for tsig object S is a row vector M whose element M(k) is
%    the standard deviation of the S's kth channel, S(k).
%  
%    Std(S,1) is the same as std(S). Thus by convention, the "first 
%    dimension" of tsig objects is the time dimension. This is consistent 
%    with tsig/vertcat appending tsig objects in time.
%
%    T=std(S,2) returns a tsig object containing the sample-wise std
%    across channels. Thus by convention, the "second dimension" of tsig
%    objects is the channel number, which is consistent with tsig/horzcat
%    concatenating channels. Importantly, the channels are time-aligned
%    before being processed. The mean tsig T runs from min([S.onset]) to 
%    max(S.offset), and for channels that do not cover this full range, 
%    zeros are padded in the missing parts.
% 
%    M=std(S,2,'nozeros') is the same as std(S,2) except that no zeros are
%    padded. In a given time range, the std will be evaluated across only 
%    those channels that exist in this range. See tsig/mean and tsig/median
%    for more details.
%
%    See also tsig/mean, tsig/meadian, tsig/plus, tsig/exist, tsig/isvoid.

if nargin<2, dim=1; end
if nargin<3, flag=''; end

if ~isequal(1,dim) && ~isequal(2,dim),
    error('Dim argument must be either 1 or 2.');
end
if isvoid(T),
    error('Invalid input: void tsig object.'); 
end
if ~isempty(flag),
    [flag, Mess] = keywordMatch(flag, {'nozeros'}, 'flag');
    error(Mess);
end

if dim==1, % std of each channel
    Y = wavefun(@std, T, 'cat');
elseif dim==2, % std across channels: delegate to xchanop
    % 'nozeros' flag determines the default value passed to xchanop
    if isequal('nozeros',flag),
        defS = 'existingchannelsonly'; % exclude non-existing channel intervals
    else,
        defS = 0; % substitute zeros where channels don't exist
    end
    my_std = @(x)std(x,0,2); % pass 2nd arg = FLAG and 3rd arg = DIM.
    [Y, Mess] = xchanop(T, my_std, defS);
    error(Mess);
end



