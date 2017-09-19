function P = unwrap(P, varargin);
% tsig/unwrap - UNWRAP for tsig objects.
% 
% See also tsig.
% 
for ii=1:nchan(P),
    P.Waveform{ii} = feval(@unwrap, P.Waveform{ii}, varargin{:});
end;
 
