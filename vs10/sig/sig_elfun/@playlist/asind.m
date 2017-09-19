function P = asind(P, varargin);
% playlist/asind - ASIND for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@asind, P.Waveform{ii}, varargin{:});
end;
 
