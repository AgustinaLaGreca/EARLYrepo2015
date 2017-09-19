function P = floor(P, varargin);
% playlist/floor - FLOOR for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@floor, P.Waveform{ii}, varargin{:});
end;
 
