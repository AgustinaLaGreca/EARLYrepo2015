function P = log2(P, varargin);
% playlist/log2 - LOG2 for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@log2, P.Waveform{ii}, varargin{:});
end;
 
