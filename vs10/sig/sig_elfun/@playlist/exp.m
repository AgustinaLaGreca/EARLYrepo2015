function P = exp(P, varargin);
% playlist/exp - EXP for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@exp, P.Waveform{ii}, varargin{:});
end;
 
