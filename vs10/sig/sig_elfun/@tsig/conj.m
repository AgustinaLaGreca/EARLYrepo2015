function P = conj(P, varargin);
% tsig/conj - CONJ for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@conj, P.Waveform{ii}, varargin{:});
end;
 
