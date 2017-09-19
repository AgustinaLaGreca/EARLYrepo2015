function P = hypot(P, varargin);
% tsig/hypot - HYPOT for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@hypot, P.Waveform{ii}, varargin{:});
end;
 
