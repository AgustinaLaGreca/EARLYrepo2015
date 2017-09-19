function P = acotd(P, varargin);
% playlist/acotd - ACOTD for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@acotd, P.Waveform{ii}, varargin{:});
end;
 
