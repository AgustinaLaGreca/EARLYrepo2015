function P = log1p(P, varargin);
% playlist/log1p - LOG1P for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@log1p, P.Waveform{ii}, varargin{:});
end;
 
