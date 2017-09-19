function P = log2(P, varargin);
% tsig/log2 - LOG2 for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@log2, P.Waveform{ii}, varargin{:});
end;
 
