function P = cplxpair(P, varargin);
% tsig/cplxpair - CPLXPAIR for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@cplxpair, P.Waveform{ii}, varargin{:});
end;
 
