function P = cosh(P, varargin);
% playlist/cosh - COSH for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@cosh, P.Waveform{ii}, varargin{:});
end;
 
