function P = acot(P, varargin);
% tsig/acot - ACOT for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@acot, P.Waveform{ii}, varargin{:});
end;
 
