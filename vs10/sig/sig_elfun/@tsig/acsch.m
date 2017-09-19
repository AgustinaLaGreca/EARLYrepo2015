function P = acsch(P, varargin);
% tsig/acsch - ACSCH for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@acsch, P.Waveform{ii}, varargin{:});
end;
 
