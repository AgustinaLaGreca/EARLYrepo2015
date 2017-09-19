function P = acsc(P, varargin);
% tsig/acsc - ACSC for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@acsc, P.Waveform{ii}, varargin{:});
end;
 
