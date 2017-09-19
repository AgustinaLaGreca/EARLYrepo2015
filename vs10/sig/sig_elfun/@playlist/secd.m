function P = secd(P, varargin);
% playlist/secd - SECD for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@secd, P.Waveform{ii}, varargin{:});
end;
 
