function P = sech(P, varargin);
% playlist/sech - SECH for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@sech, P.Waveform{ii}, varargin{:});
end;
 
