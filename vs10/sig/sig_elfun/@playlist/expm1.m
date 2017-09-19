function P = expm1(P, varargin);
% playlist/expm1 - EXPM1 for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@expm1, P.Waveform{ii}, varargin{:});
end;
 
