function P = cscd(P, varargin);
% playlist/cscd - CSCD for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@cscd, P.Waveform{ii}, varargin{:});
end;
 
