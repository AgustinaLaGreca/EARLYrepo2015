function P = acotd(P, varargin);
% tsig/acotd - ACOTD for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@acotd, P.Waveform{ii}, varargin{:});
end;
 
