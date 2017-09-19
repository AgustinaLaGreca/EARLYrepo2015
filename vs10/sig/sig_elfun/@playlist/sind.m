function P = sind(P, varargin);
% playlist/sind - SIND for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@sind, P.Waveform{ii}, varargin{:});
end;
 
