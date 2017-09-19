function P = asech(P, varargin);
% tsig/asech - ASECH for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@asech, P.Waveform{ii}, varargin{:});
end;
 
