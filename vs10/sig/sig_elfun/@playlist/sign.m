function P = sign(P, varargin);
% playlist/sign - SIGN for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@sign, P.Waveform{ii}, varargin{:});
end;
 
