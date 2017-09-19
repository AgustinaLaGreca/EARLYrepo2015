function [Y, varargout] = wavefun(Fun, T, flag, varargin);
% tsig/wavefun - apply a function to the waveforms of a tsig object
%   C=wavefun(FUN,T) applies function FUN to the waveforms of tsig T and
%   returns the result in cell array C. The cells of C correspond to the
%   channels of T:  C{ichan} = FUN(T.Waveform{k}). FUN is a function handle.
%
%   [C,D, ..]=wavefun(FUN,T) returns any additional outputs provided by FUN.
%
%   C=wavefun(FUN,T,'cell') is the same as C=wavefun(FUN,T).
%
%   [S,D,..]=wavefun(FUN,T,'tsig') returns a tsig object containing the results of
%   the channelwise FUN calls as waveforms. By convention the sample rate
%   and onsets of S are identical to those of T.
%
%   [C,D,..]=wavefun(FUN,T,'cat') concatenates the cells C{ichan},D{ichan} ..
%   of the different channels. This call is equivalent to
%         [C,D,..]=wavefun(FUN,T,'cell'); % outputs are cell arrays C,D..
%         C = [C{:}]; % paste cells into single array.
%         D = [D{:}];
%         ...
%   Of course, the outputs C,D,.. must be uniform for this to succeed.
%
%   C=wavefun(FUN,T,'cell', arg1, ..) and C=wavefun(FUN,T,'cat', arg1, ..)
%   pass the additional arguments arg1, .. to the FUN calls. That is, 
%   FUN is called as  FUN(T.Waveform{k}, arg1, ..)
%
%   Examples
%     C = wavefun(@max, T) returns the channelwise maxima in cell array C
%     M = wavefun(@max, T, 'cat') returns the maxima in row vector M
%     T = wavefun(@max, T, 'tsig', 0) rectifies the waveforms of T.
%
%   See also tsig/subsref, tsig/find, CELLUN.

if ~isa(Fun, 'function_handle'),
    error('First argument is not a function handle.');
end
if nargin<3, flag='cell'; end

[flag, Mess] = keywordMatch(flag, {'cell' 'cat' 'tsig'}, 'flag');
error(Mess);
% convert flag into cellfun modifier
UniformOutput = isequal('cat', flag);

varargout = cell(1,nargout-1); % initialize so cellfun knows how many outargs to provide
fun = @(x)Fun(x,varargin{:}); % anonymous fnc to automatically incorporate additional args
[Y, varargout{:}] = cellfun(fun, T.Waveform, 'UniformOutput', UniformOutput);

if isequal('tsig', flag), % re-convert top tsig. Use constructor to intercept any illegal contents
    Y = tsig(fsam(T), Y, onset(T));
end




