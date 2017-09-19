function P = complex(P, varargin);
% playlist/complex - COMPLEX for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@complex, P.Waveform{ii}, varargin{:});
end;
 
