function P = sqrt(P, varargin);
% playlist/sqrt - SQRT for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@sqrt, P.Waveform{ii}, varargin{:});
end;
 
