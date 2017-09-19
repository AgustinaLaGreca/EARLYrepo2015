function P = asinh(P, varargin);
% playlist/asinh - ASINH for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@asinh, P.Waveform{ii}, varargin{:});
end;
 
