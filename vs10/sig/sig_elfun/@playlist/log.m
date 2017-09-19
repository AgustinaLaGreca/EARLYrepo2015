function P = log(P, varargin);
% playlist/log - LOG for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@log, P.Waveform{ii}, varargin{:});
end;
 
