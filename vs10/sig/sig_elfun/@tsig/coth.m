function P = coth(P, varargin);
% tsig/coth - COTH for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@coth, P.Waveform{ii}, varargin{:});
end;
 
