function P = asecd(P, varargin);
% playlist/asecd - ASECD for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@asecd, P.Waveform{ii}, varargin{:});
end;
 
