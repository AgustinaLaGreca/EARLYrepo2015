function P = log1p(P, varargin);
% tsig/log1p - LOG1P for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@log1p, P.Waveform{ii}, varargin{:});
end;
 
