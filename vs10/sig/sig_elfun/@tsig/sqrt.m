function P = sqrt(P, varargin);
% tsig/sqrt - SQRT for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@sqrt, P.Waveform{ii}, varargin{:});
end;
 
