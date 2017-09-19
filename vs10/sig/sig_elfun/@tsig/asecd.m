function P = asecd(P, varargin);
% tsig/asecd - ASECD for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@asecd, P.Waveform{ii}, varargin{:});
end;
 
