function P = cplxpair(P, varargin);
% playlist/cplxpair - CPLXPAIR for playlist objects.
% 
% See also playlist.
% 
for ii=1:nwave(P),
    P.Waveform{ii} = feval(@cplxpair, P.Waveform{ii}, varargin{:});
end;
 
