function P = sech(P, varargin);
% tsig/sech - SECH for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@sech, P.Waveform{ii}, varargin{:});
end;
 
