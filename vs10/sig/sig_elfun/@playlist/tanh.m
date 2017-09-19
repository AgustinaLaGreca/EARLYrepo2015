function P = tanh(P, varargin);
% playlist/tanh - TANH for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@tanh, P.Waveform{ii}, varargin{:});
end;
 
