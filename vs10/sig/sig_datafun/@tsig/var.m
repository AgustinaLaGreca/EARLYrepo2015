function Y = var(T, dim,flag);
% tsig/var - VAR for tsig objects.
%    M=var(S) for tsig object S is a row vector M whose element M(k) is
%    the variance of the S's kth channel, S(k).
%  
%    var(S,1) is the same as var(S). Thus by convention, the "first 
%    dimension" of tsig objects is the time dimension. This is consistent 
%    with tsig/vertcat appending tsig objects in time.
%
%    T=var(S,2) returns a tsig object containing the sample-wise variance
%    across channels. Thus by convention, the "second dimension" of tsig
%    objects is the channel number, which is consistent with tsig/horzcat
%    concatenating channels. Importantly, the channels are time-aligned
%    before being processed. The tsig T runs from min([S.onset]) to 
%    max(S.offset), and for channels that do not cover this full range, 
%    zeros are padded in the missing parts.
% 
%    M=var(S,2,'nozeros') is the same as var(S,2) except that no zeros are
%    padded. In a given time range, the var will be evaluated across only 
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

if dim==1, % var of each channel
    Y = wavefun(@var, T, 'cat');
elseif dim==2, % var across channels: delegate to xchanop
    % 'nozeros' flag determines the default value passed to xchanop
    if isequal('nozeros',flag),
        defS = 'existingchannelsonly'; % exclude non-existing channel intervals
    else,
        defS = 0; % substitute zeros where channels don't exist
    end
    my_var = @(x)var(x,0,2); % pass 2nd arg = W and 3rd arg = DIM.
    [Y, Mess] = xchanop(T, my_var, defS);
    error(Mess);
end



