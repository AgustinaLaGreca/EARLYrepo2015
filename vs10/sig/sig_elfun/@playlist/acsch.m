function P = acsch(P, varargin);
% playlist/acsch - ACSCH for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@acsch, P.Waveform{ii}, varargin{:});
end;
 
