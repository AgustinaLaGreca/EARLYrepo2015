function P = log10(P, varargin);
% playlist/log10 - LOG10 for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@log10, P.Waveform{ii}, varargin{:});
end;
 
