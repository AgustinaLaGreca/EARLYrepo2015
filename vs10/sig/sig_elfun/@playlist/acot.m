function P = acot(P, varargin);
% playlist/acot - ACOT for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@acot, P.Waveform{ii}, varargin{:});
end;
 
