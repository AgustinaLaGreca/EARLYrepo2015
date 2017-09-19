function P = reallog(P, varargin);
% playlist/reallog - REALLOG for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@reallog, P.Waveform{ii}, varargin{:});
end;
 
