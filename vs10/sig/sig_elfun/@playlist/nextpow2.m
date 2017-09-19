function P = nextpow2(P, varargin);
% playlist/nextpow2 - NEXTPOW2 for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@nextpow2, P.Waveform{ii}, varargin{:});
end;
 
