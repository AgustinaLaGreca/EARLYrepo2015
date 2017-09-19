function P = pow2(P, varargin);
% playlist/pow2 - POW2 for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@pow2, P.Waveform{ii}, varargin{:});
end;
 
