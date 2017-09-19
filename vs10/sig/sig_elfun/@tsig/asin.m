function P = asin(P, varargin);
% tsig/asin - ASIN for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@asin, P.Waveform{ii}, varargin{:});
end;
 
