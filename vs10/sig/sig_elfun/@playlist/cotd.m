function P = cotd(P, varargin);
% playlist/cotd - COTD for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@cotd, P.Waveform{ii}, varargin{:});
end;
 
