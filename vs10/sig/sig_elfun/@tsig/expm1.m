function P = expm1(P, varargin);
% tsig/expm1 - EXPM1 for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@expm1, P.Waveform{ii}, varargin{:});
end;
 
