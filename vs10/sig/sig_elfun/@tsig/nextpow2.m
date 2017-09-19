function P = nextpow2(P, varargin);
% tsig/nextpow2 - NEXTPOW2 for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@nextpow2, P.Waveform{ii}, varargin{:});
end;
 
