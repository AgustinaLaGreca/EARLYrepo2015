function P = acscd(P, varargin);
% tsig/acscd - ACSCD for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@acscd, P.Waveform{ii}, varargin{:});
end;
 
