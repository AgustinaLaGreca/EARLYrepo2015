function P = angle(P, varargin);
% playlist/angle - ANGLE for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@angle, P.Waveform{ii}, varargin{:});
end;
 
