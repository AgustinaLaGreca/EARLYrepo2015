function Y = vertcat(varargin)
% tsig/vertcat - append tsig objects in time.
%   vertcat(S,T,..) or [S;T;..] is a tsig object obtained by appending, in
%   orde, S,T, .. (S then T, etc). Any onsets of the input signals are 
%   ignored. S, T .. must have compatible channel counts and identical 
%   sample frequency.
%
%   By definition, [[];S] == [S;[]] == S. Also [S; E]== [E; S] == S when E
%   is a void tsig.
%
%   See also tsig/vertcat, tsig/plus, tsig/fliplr, tsig/isvoid.

if nargin<1, % return void tsig
    Y = tsig; 
    return;
end

if ~all([cellfun(@isTsig, varargin) | cellfun(@isEmptynum, varargin)]),
    error('Tsig objects can only be appended to other Tsig objects or to [].')
end

W = {[]};
Fs = [];
for ii=1:nargin,
    if isTsig(varargin{ii}) && ~isvoid(varargin{ii}),
       Fs = [Fs fsam(varargin{ii})];
       Fs = Uniquify(Fs);
       if ~isscalar(Fs),
           error('Tsig objects can only be appended if they have the same sample frequency.');
       end
       Wapp = varargin{ii}.Waveform;
       [W, Wapp] = SameSize(W,Wapp);
       for ichan=1:length(W),
           W{ichan} = [W{ichan}; Wapp{ichan}];
       end
    end
end
Y = tsig(Fs, W); 




