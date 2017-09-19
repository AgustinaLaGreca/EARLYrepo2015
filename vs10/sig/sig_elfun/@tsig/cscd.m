function P = cscd(P, varargin);
% tsig/cscd - CSCD for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@cscd, P.Waveform{ii}, varargin{:});
end;
 
