function Y = horzcat(varargin)
% tsig/horzcat - concatenate channels of tsig objects
%   horzcat(S,T,..) or [S,T,..] is a tsig object containing the channels of 
%   S, T, .. . All arguments must be tsig objects, with the exception of
%   []. By convention [S []] == [[] S] == S. Likewise, [S E] = S when E is
%   a void tsig. All tsigs must have the same smaple frequency.
%
%   See tsig/vertcat, tsig/fliplr, tsig/isvoid.

Y = tsig;

for ii=1:nargin,
    if isempty(varargin{ii}), continue; end
    if ~isTsig(varargin{ii}),
        error('Tsig objects can only be concatenated with other tsig objects or [].');
    end
    Y.Fsam = [Y.Fsam, varargin{ii}.Fsam]; 
    Y.Waveform = [Y.Waveform, varargin{ii}.Waveform]; 
    Y.t0 = [Y.t0, varargin{ii}.t0]; 
end
Y.Fsam = Uniquify(Y.Fsam);
if ~isscalar(Y.Fsam),
    error('Tsig objects having different sample frequencies cannot be concatenated.');
end




