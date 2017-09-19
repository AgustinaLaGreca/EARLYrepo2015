function P = complex(P, varargin);
% tsig/complex - COMPLEX for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@complex, P.Waveform{ii}, varargin{:});
end;
 
