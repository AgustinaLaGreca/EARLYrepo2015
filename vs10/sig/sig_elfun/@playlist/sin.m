function P = sin(P, varargin);
% playlist/sin - SIN for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@sin, P.Waveform{ii}, varargin{:});
end;
 
