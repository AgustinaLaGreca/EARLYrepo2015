function P = asech(P, varargin);
% playlist/asech - ASECH for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@asech, P.Waveform{ii}, varargin{:});
end;
 
