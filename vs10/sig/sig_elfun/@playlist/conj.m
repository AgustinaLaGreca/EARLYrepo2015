function P = conj(P, varargin);
% playlist/conj - CONJ for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@conj, P.Waveform{ii}, varargin{:});
end;
 
