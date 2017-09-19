function P = atand(P, varargin);
% playlist/atand - ATAND for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@atand, P.Waveform{ii}, varargin{:});
end;
 
