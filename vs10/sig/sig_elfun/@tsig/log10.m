function P = log10(P, varargin);
% tsig/log10 - LOG10 for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@log10, P.Waveform{ii}, varargin{:});
end;
 
