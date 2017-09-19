function P = mod(P, varargin);
% playlist/mod - MOD for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@mod, P.Waveform{ii}, varargin{:});
end;
 
