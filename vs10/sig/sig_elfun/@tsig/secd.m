function P = secd(P, varargin);
% tsig/secd - SECD for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@secd, P.Waveform{ii}, varargin{:});
end;
 
