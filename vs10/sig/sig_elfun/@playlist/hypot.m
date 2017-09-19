function P = hypot(P, varargin);
% playlist/hypot - HYPOT for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@hypot, P.Waveform{ii}, varargin{:});
end;
 
